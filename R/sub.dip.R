sub.dip <- function(N.ITERS,VEC.IN,SUB.SAMPLE.SIZE) {
    #auxiliary function used for sub-sampling the dip statistic.
    #N.ITERS is the number of sub-sampling iterations
    #VEC.IN is the "population" of data from which we will compute the dip statistic.
    #SUB.SAMPLE.SIZE is the number of observations we will sample (with replacement) from VEC.IN 
    #the function dip.test is form the diptest library
    require(diptest)
    emp <- rep(NA,N.ITERS)
    for (i in seq(N.ITERS)) {
        NEGATIVES <- TRUE
        while (NEGATIVES) {
            sub.s <- sample(VEC.IN,size=SUB.SAMPLE.SIZE,replace=TRUE)
            if (min(sub.s) > 0) NEGATIVES <- FALSE
        }
        emp[i] <- suppressMessages(dip.test(sub.s))$p.value
    }
    return(mean(emp))
}
