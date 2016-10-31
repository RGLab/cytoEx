#' determine the best channel according to a two-stage process, described below
#'
#' STEP 1: is the same as best.dip.R
#'
#' STEP 2: If a unique minimum p-value exists, select the corresponding channel.
#'         Otherwise, sub-sample from each channel and create distribution of p-values for each candidate channel (they tied).
#'         Compute mean of each sub-sample.
#'         Initially, consider only channels with sub-sampled means below threshold := 0.05
#'         If no channels meet this criteria, increment the threshold by 0.01 until at least one channel is in contention.
#'         If a unique channel falls below the threshold, pick it.
#'         If multiple channels fall below the threshold, pick model with maximum difference in ICL, where the difference is:
#'            (ICL from mixture with two components) - (ICL from mixture with three components).
#' @importFrom diptest dip.test
#' @importFrom flowClust flowClust criterion
#' @importFrom parallel mclapply
#' @param ALPHA user selected significance level for channel-wide dip statistic. bonferonni correction applied to this value.
#' @param P.ITERS number of sub-samples taken from channels which pass initial screen
#' @param SS.SIZE size of sub-sample taken from channel which passes the initial screen
#' @return Data table with decision metrics for gating method.
best.dip.ICL <- function(fr, debug.mode=FALSE, plotEnv=new.env(parent=emptyenv()), parallel_type, mc.cores,
                         ALPHA = 0.05, P.ITERS=10000, SS.SIZE = 200, ...) {

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
    #first.screen <- potential.channels[intersect(which(first.p.list == min(unique(first.p.list))),which(first.p.list < bonferonni.alpha))]
    first.screen <- potential.channels[which(first.p.list < bonferonni.alpha)]
    first.remainders <- setdiff(potential.channels,first.screen)

    #record ranks for reporting/visualization
    #scores <- rep((length(potential.channels)+1),length(potential.channels))
    scores <- rep(0,length(potential.channels))
    CURRENT.CHANNEL.RANK <- 1
    for (chn in potential.channels) {
        if (chn %in% first.remainders) {
            scores[which(potential.channels==chn)] <- CURRENT.CHANNEL.RANK
            CURRENT.CHANNEL.RANK <- CURRENT.CHANNEL.RANK + 1
        }
    }

    # if no channels pass the intial screening, return a table with  all scores set to -1
    if (length(first.screen) == 0) {
        return.table <- data.table(channel=potential.channels,score=-1,area.ratio=0,intial.p.values = first.p.list, second.p.values = 1, diff.icl = -Inf, b.alpha = bonferonni.alpha)
        return(return.table)
    }

    # if only one channel passes the initial screening, pick it.
    else  if (length(first.screen) == 1){
        second.pv <- rep(1, length(potential.channels)) #this vector is for reporting
        second.pv[which(potential.channels==first.screen)] <- min(unique(first.p.list))
        selected.channel <- first.screen
        icl.record <- rep(-Inf, length(potential.channels)) #this vector is also for reporting
        mixtures <- flowClust(fr,varNames=c(first.screen),K=2:3, B=10000, lambda=1, trans=0)
        plotObjs <- list(ggflowClust.hist(mixtures[[1]], fr))
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

            mixtures <- flowClust(fr,varNames=c(CAND),K=2:3
                                  # , B=10000
                                  , lambda=1
                                  , trans = 0, min.count = -1, max.count = -1, randomStart = 0)
            m.icls <- criterion(mixtures,"ICL")
            icl.difference <- m.icls[1]-m.icls[2]
            if(debug.mode){
                p <- ggflowClust.hist(mixtures[[1]], fr)
            }else
                p <- NULL

            if (parallel_type == "none") dip.sub <- mean(unlist(lapply(emp,FUN=sub.sample.dip.p)))
            else dip.sub <- mean(unlist(parallel::mclapply(emp,FUN=sub.sample.dip.p)))
            return(list(scores=list(dip.sub,icl.difference),plot=p))
        }

        #joint list is a flattened list matched to candidate channel in the following way
        #odd index: the mean of a sub-sampled p.value for a given channel
        #odd index + 1: the difference in icl for that channel
        if (parallel_type == "none") res <- lapply(first.screen,FUN=get.dip.and.icl)
        else res <- parallel::mclapply(first.screen,FUN=get.dip.and.icl, mc.cores = mc.cores)
        joint.list <- unlist(sapply(res, function (x) {unlist(x[["scores"]])}))
        plotObjs <- sapply(res, function (x) {x[["plot"]]}, simplify = FALSE)

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

        #finally, use sub-sampled p-values determine which channels remain in contention.
        thresh.p <- 0.05
        NULL.SCREEN <- TRUE
        while (NULL.SCREEN) {
            second.screen <- first.screen[which(p.list <= thresh.p)]
            if (length(second.screen) == 0) {
                thresh.p <- thresh.p + 0.01
            }
            else {
                NULL.SCREEN <- FALSE
            }
        }

        #update rank vector
        second.remainders <- setdiff(first.screen,second.screen)
        for (chn in potential.channels) {
            if (chn %in% second.remainders) {
                scores[which(potential.channels==chn)] <- CURRENT.CHANNEL.RANK
                CURRENT.CHANNEL.RANK <- CURRENT.CHANNEL.RANK + 1
            }
        }

        #get remaining channels difference of icls
        dicl.sub <- dicl.list[which(p.list <= thresh.p)]

        #pick channel with largest difference (indicating we prefer it)
        selected.channel <- second.screen[which(dicl.sub == max(unique(dicl.sub)))]

        #final ranks
        repdf <- data.frame(channels=second.screen,scores=dicl.sub)
                                        #repdf <- repdf[order(repdf$scores,decreasing=TRUE),]
        repdf <- repdf[order(repdf$scores),]
        for (chn in repdf$channels) {
            scores[which(potential.channels==chn)] <- CURRENT.CHANNEL.RANK
            CURRENT.CHANNEL.RANK <- CURRENT.CHANNEL.RANK + 1
        }

    }

    #at this point, do not expect numerical ties for p.values because we sub-sampled.
    if (length(selected.channel) > 1) {
        stop("More than one channel selected in best.dip call")
    }

    if (debug.mode) {
        j <- 1
        for(cn in potential.channels){
            if (!(cn %in% first.screen)) {
                mixtures <- flowClust(fr,varNames=c(cn),K=2:3, B=10000, lambda=1, trans=0)
                p <- ggflowClust.hist(mixtures[[1]], fr) + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
            } else {
                p <- plotObjs[[j]] + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
                j <- j+1
            }
            p <- p + scale_colour_discrete(guide = FALSE) + ylab("")
            marker <- getChannelMarker(fr, cn)[["desc"]]
            p <- p + xlab(marker)
            plotEnv[[cn]] <- p
        }
    }

    #create a vector for reporting compatibility with cytoEx
    scores[which(potential.channels==selected.channel)] <- length(potential.channels)
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
