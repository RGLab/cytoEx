#' @importFrom diptest dip.test
#' @importFrom parallel mclapply
#'
#' determine the best channel according to a two-stage process, described below
#'        
#' STEP 1: for each channel, compute p-values for Hartigan's dip statistic (with bonferonni correction) to decide if any channels are multimodal. 
#'        Threshold at \alpha = ALPHA/(# of channels) 
#'        If no channel is below \alpha, return data.table with all areas == 0 (indicates no channel selected to ambient gating function)
#'        Else, proceed to step 2.
#'
#' STEP 2: If a unique minimum p-value exists, select the corresponding channel.
#'        Otherwise, sub-sample from each channel and create distribution of dip-statistic p-values for each candidate channel (those that tied). 
#'        Compute mean of each sub-sample, and select channel with minimum average p-value.
#'
#' @param ALPHA user selected significance level for channel-wide dip statistic. bonferonni correction applied to this value.
#' @param P.ITERS number of sub-samples taken from channels which pass initial screen
#' @param SS.SIZE size of sub-sample taken from channel which passes the initial screen
#'
#' @return Data table with decision metrics for gating method.
best.dip <- function(fr, debug.mode=FALSE, plotEnv=new.env(parent=emptyenv()), parallel_type, mc.cores,
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
        first.p.list <- append(first.p.list,suppressMessages(diptest::dip.test(cyto.data[,which(colnames(cyto.data) == candidate)]))$p.value)
    }
    
    #channels which pass the first screen have p values below the bonferonni-adjusted significance level.
    first.screen <- potential.channels[intersect(which(first.p.list == min(unique(first.p.list))),which(first.p.list < bonferonni.alpha))]
    
    # if no channels pass the intial screening, return a table with  all scores set to -1
    if (length(first.screen) == 0) {
        return.table <- data.table(channel=potential.channels,score=-1,area.ratio=0,intial.p.values = first.p.list, second.p.values = 1, b.alpha = bonferonni.alpha)
        return(return.table)
    }

    # if only one channel passes the initial screening, pick it.
    else  if (length(first.screen) == 1){
        second.pv <- rep(1, length(potential.channels)) #this vector is for reporting
        second.pv[which(potential.channels==first.screen)] <- min(unique(first.p.list))
        selected.channel <- first.screen
    }

    # finally, if more than one channel passes the intial screen, there are ties. proceed to stage two.
    else {
        second.pv <- rep(1, length(potential.channels)) #this vector is for reporting

        get.sub.sample.dip.p <- function(CAND) {
            cand.data <- cyto.data[,which(colnames(cyto.data) == CAND)]
            emp <- rep(NA,P.ITERS)
            sub.sample.dip.p <- function(x) {
                NEGATIVES <- TRUE
                while (NEGATIVES) {
                    sub.s <- sample(cand.data,size=SS.SIZE,replace=FALSE)
                    if (min(sub.s) > 0) NEGATIVES <- FALSE
                }
                return(suppressMessages(diptest::dip.test(sub.s))$p.value)
            }
            if (parallel_type == "none") return(mean(unlist(lapply(emp,FUN=sub.sample.dip.p))))
            else return(mean(unlist(parallel::mclapply(emp,FUN=sub.sample.dip.p))))
        }
        
        if(parallel_type == "none") p.list <- unlist(lapply(first.screen,FUN=get.sub.sample.dip.p))
        else p.list <- unlist(parallel::mclapply(first.screen,FUN=get.sub.sample.dip.p, mc.cores = mc.cores))
                                            
        #record the mean of sub-sampled p-values
        i <- 1
        for (candidate in first.screen) {
            second.pv[which(potential.channels==candidate)] <- p.list[i]
            i <- i+1
        }
        selected.channel <- first.screen[which(p.list == min(unique(p.list)))]
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
    return.table <- data.table(channel=potential.channels,score=scores,area.ratio=areas,intial.p.values = first.p.list, second.p.values=second.pv,b.alpha = bonferonni.alpha)
    return(return.table)
}
