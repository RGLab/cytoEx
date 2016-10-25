#' @importFrom flowClust flowClust
best.dip.ICL <- function(fr, debug.mode=FALSE, plotEnv=new.env(parent=emptyenv()), parallel_type, mc.cores,
                     CHANNELS = c(), ALPHA = 0.01, P.ITERS=10000, SS.SIZE = 200, ...) {
    #determine the best channel according to a two-stage process, described below
        
    #STEP 1: is the same as best.dip.R

    #STEP 2: If a unique minimum p-value exists, select the corresponding channel.
    #        Otherwise, sub-sample from each channel and create bootstrap distribution of p-values for each candidate channel (those with same minimum). 
    #        Compute mean of each bootstrap sample.
    #        Initially, consider only channels with bootstrap means below threshold 0.05 
    #        If no channels meet this criteria, increment the threshold by 0.01 until at least one channel is in contention.
    #        If a unique channel falls below the threshold, pick it.
    #        If multiple channels fall below the threshold, pick model with maximum difference in ICL, where the difference is:
    #           (ICL from mixture with two components) - (ICL from mixture with three components).
    require(flowClust)
    require(diptest)
    #if the caller provides a customized vector of potential channels, use it.
    if (length(CHANNELS) > 0) {
        potential.channels <-  CHANNELS
    }
    #otherwise, grab all the channels from the flow frame, and then remove forward scatter, side scatter, and time columns (if present).
    else {
        bad.channels <- c("FSC-A","FSC-H","FSC-W","SSC-A","Time")
        all.channels <- fr@parameters$name
        potential.channels <- setdiff(all.channels,bad.channels)
    }
    #with the potential channels selected, apply bonferonni correction
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
        icl.record <- second.pv <- rep(-Inf, length(potential.channels)) #this vector is for reporting
        dicl.list <- p.list <- c()
        for (candidate in first.screen) {
            sub.sampled.p.value <- sub.dip(P.ITERS,cyto.data[,which(colnames(cyto.data) == candidate)],SS.SIZE)
            second.pv[which(potential.channels==candidate)] <- sub.sampled.p.value
            p.list <- append(p.list,sub.sampled.p.value)
            mixtures <- flowClust(fr,varNames=c(candidate),K=2:3, B=10000, lambda=1, trans=0)
            m.icls <- criterion(mixtures,"ICL")
            icl.difference <- m.icls[1]-m.icls[2]
            icl.record[which(potential.channels==candidate)] <- icl.difference
            dicl.list <- append(dicl.list,icl.difference)
        }
        
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
        dicl.sub <- dicl.list[which(p.list <= thresh.p)]
        selected.channel <- second.screen[which(dicl.sub == max(unique(dicl.sub)))]
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
