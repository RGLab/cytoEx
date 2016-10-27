#' determine the best channel according to a two-stage process, described below
#'    
#' STEP 1: is the same as best.dip.R
#'
#' STEP 2: If a unique minimum p-value exists, select the corresponding channel.
#'         Otherwise, sub-sample from each channel and create a distribution of p-values for each candidate channel (they tied).
#'         Compute mean of each bootstrap sample AND difference in ICL for each channel where, for a given channel, we compute
#'            (ICL from mixture with two components) - (ICL from mixture with three components).
#'         Subset to only consider models with a positive difference in ICL
#'         Among these models, pick the model with the minimum mean sub-sampled p-value.
#' @importFrom diptest dip.test
#' @importFrom flowClust flowClust
#' @importFrom parallel mclapply 
#' @param ALPHA user selected significance level for channel-wide dip statistic. bonferonni correction applied to this value.
#' @param P.ITERS number of sub-samples taken from channels which pass initial screen
#' @param SS.SIZE size of sub-sample taken from channel which passes the initial screen
#' @return Data table with decision metrics for gating method.

best.dip.ICL.sub <- function(fr, debug.mode=FALSE, plotEnv=new.env(parent=emptyenv()), parallel_type, mc.cores,
                             ALPHA = 0.01, P.ITERS=10000, SS.SIZE = 200, ...) {
    potential.channels <- fr@parameters$name
    #apply bonferonni correction
    bonferonni.alpha <- ALPHA/length(potential.channels)

    #get the cytometry data from the flow frame
    cyto.data <- exprs(fr)

    #conduct stage one of the channel selection procedure
    first.p.list <- c()
    for (candidate in potential.channels) {
        #suppress warnings because can have greater than 70000 observations
        first.p.list <- append(first.p.list,suppressMessages(dip.test(cyto.data[,which(colnames(cyto.data) == candidate)]))$p.value)
    }
    
    #channels which pass the first screen have p values below the bonferonni-adjusted significance level.
    first.screen <- potential.channels[intersect(which(first.p.list == min(unique(first.p.list))),which(first.p.list < bonferonni.alpha))]
    
    # if no channels pass the intial screening, return a table with  all scores set to -1
    if (length(first.screen) == 0) {
        return.table <- data.table(channel=potential.channels,score=-1,area.ratio=0,intial.p.values = first.p.list, second.p.values = 1, icl.record = -Inf, b.alpha = bonferonni.alpha)
        return(return.table)
    }

    # if only one channel passes the initial screening, pick it.
    else  if (length(first.screen) == 1){
        second.pv <- rep(1, length(potential.channels)) #this vector is for reporting
        second.pv[which(potential.channels==first.screen)] <- min(unique(first.p.list))
        selected.channel <- first.screen
        icl.record <- rep(-Inf, length(potential.channels)) #this vector is also for reporting
    }

    # finally, if more than one channel passes the intial screen, there are ties. proceed to stage two.
    else {
        get.dip.and.icl <- function(CAND) {
            cand.data <- cyto.data[,which(colnames(cyto.data) == CAND)]
            emp <- rep(NA,P.ITERS)
            
            sub.sample.dip.p <- function(x) {
                NEGATIVES <- TRUE
                while (NEGATIVES) {
                    sub.s <- sample(cand.data,size=SS.SIZE,replace=TRUE)
                    if (min(sub.s) > 0) NEGATIVES <- FALSE
                }
                return(suppressMessages(diptest::dip.test(sub.s))$p.value)
            }
            
            mixtures <- flowClust(fr,varNames=c(CAND),K=2:3, B=10000, lambda=1, trans=0)
            m.icls <- criterion(mixtures,"ICL")
            icl.difference <- m.icls[1]-m.icls[2]
            
            if (parallel_type == "none") dip.sub <- mean(unlist(lapply(emp,FUN=sub.sample.dip.p)))
            else dip.sub <- mean(unlist(parallel::mclapply(emp,FUN=sub.sample.dip.p)))
            return(list(dip.sub,icl.difference))
        }

        #joint list is a flattened list matched to candidate channel in the following way
        #odd index: the mean of a sub-sampled p.value for a given channel
        #odd index + 1: the difference in icl for that channel
        if (parallel_type == "none") joint.list <- unlist(lapply(first.screen,FUN=get.dip.and.icl))
        else joint.list <- unlist(parallel::mclapply(first.screen,FUN=get.dip.and.icl, mc.cores = mc.cores))

        #record these quanities for metrics, and generate p.list and  mean of sub-sampled p-values
        second.pv <- rep(1, length(potential.channels)) 
        icl.record <- rep(-Inf, length(potential.channels))
        dicl.list <- p.list <- c()
        i <- 1
        for (candidate in first.screen) {
            second.pv[which(potential.channels==candidate)] <- joint.list[i]
            icl.record[which(potential.channels==candidate)] <- joint.list[(i+1)]
            p.list <- append(p.list,joint.list[i])
            dicl.list <- append(dicl.list,joint.list[(i+1)])
            i <- (i+2) #step by two to arrive at next sub-sampled p.value.
        }
        
        
        #first, sub-set to models that seem best described by a 2-component mixture 
        second.screen <- first.screen[which(dicl.list >= 0)]
        if (length(second.screen) == 0) {
            stop("Second screen appears to be better fit by a 3-component mixture, using difference of ICLs")
        }
        
        #then, use the sub-sampled p-value to select among remaning candidates
        p.sub <- p.list[which(dicl.list >= 0)]
        selected.channel <- second.screen[which(p.sub == min(unique(p.sub)))]
        
    }
    

    #at this point, do not expect numerical ties for p.values because we sub-sampled.
    if (length(selected.channel) > 1) {
        stop("More than one channel selected in best.dip call")
    }

    #create a vector for reporting compatibility with cytoEx
    scores <- rep(-1, length(potential.channels))
    scores[which(potential.channels==selected.channel)] <- 1
    areas <- rep(0, length(potential.channels))
    areas[which(potential.channels==selected.channel)] <- 1
    return.table <- data.table(channel=potential.channels,
                               score=scores,
                               area.ratio=areas,
                               intial.p.values = first.p.list,
                               second.p.values=second.pv,
                               diff.icl=icl.record,
                               b.alpha = bonferonni.alpha)
    return(return.table)
}
