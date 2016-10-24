sub.dip <- function(N.ITERS,VEC.IN,SUB.SAMPLE.SIZE) {
    #auxiallary function for best.dip. used to compute dip statistic on sub-samples.
    #N.ITERS is the number of sub-sampling iterations
    #VEC.IN is the "population" of data from which we will compute the dip statistic.
    #SUB.SAMPLE.SIZE is the number of observations we will sample (with replacement) from VEC.IN 
    #the function dip.test is form the diptest library
    require(diptest)
    emp <- rep(NA,N.ITERS)
    for (i in seq(N.ITERS)) {
        sub.s <- sample(VEC.IN,size=SUB.SAMPLE.SIZE,replace=TRUE)
        emp[i] <- suppressMessages(dip.test(sub.s))$p.value
    }
    return(mean(emp))
}

best.dip <- function(fr, debug.mode=FALSE, plotEnv=new.env(parent=emptyenv()), parallel_type, mc.cores,
                     CHANNELS = c(), ALPHA = 0.01, P.ITERS=10000, ...) {
    #determine the best channel according to a two-stage process, described below
        
    #STEP 1: for each channel, compute p-values for Hartigan's dip statistic (with bonferonni correction) to decide if any channels are multimodal. 
    #        Threshold at \alpha = ALPHA/(# of channels) 
    #        If no channel is below \alpha, return data.table with all areas == 0 (currently indicates no channel selected to ambient gating function)
    #        Else, proceed to step 2.

    #STEP 2: If a unique minimum p-value exists, select the corresponding channel.
    #        Otherwise, sub-sample from each channel and create bootstrap distribution of p-values for each candidate channel (those with same minimum). 
    #        Compute mean of each bootstrap, and select channel with minimum average p-value.
    #
    #        Intuition: channels which exhibit distinct modes of roughly equal mass are channels which investigators would notice first when clustering by hand.
    #                   sub-samples of such channels will have, on average, small dip-statistics since each mode will tend to appear in a sub-sample.
    #                   on the other hand, if a channel has multiple modes but one is predominant, sub-samples will tend to look unimodal (and so have larger dip statistic p-values).
    #                   method assumes that the channels with clear separation would be gated first by invesitgators.
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
        p.list <- c()
        for (candidate in first.screen) {
            sub.sampled.p.value <- sub.dip(P.ITERS,cyto.data[,which(colnames(cyto.data) == candidate)],200)
            second.pv[which(potential.channels==candidate)] <- sub.sampled.p.value
            p.list <- append(p.list,sub.sampled.p.value)
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
