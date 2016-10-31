library(flowWorkspace)


test_that("Test best.dip singlet selection unchanged: ", {
    set.seed(1239214)
    dataDir <- system.file("extdata",package="flowWorkspaceData")
    gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
    gs.dip <- clone(gs, isNew = FALSE, isEmpty = FALSE)
    for(node in getChildren(gs.dip, "singlets"))
        Rm(node, gs.dip)
    fs <- getData(gs.dip,"singlets") 
    fs <- getData(gs.dip, "singlets")
    fr <- fs[[1, use.exprs = FALSE]]
    pd <- pData(parameters(fr))
    pd <- pd[!is.na(pd[["desc"]]),]
    channels <- pd[["name"]]
    fr <- fs[[1,channels]] #we now only have fluorescent channels
    res <- cytoEx:::best.dip.ICL(fr,parallel_type="multicore",mc.cores=4)
    
    #this data set should have two channels with non-zero p-value for dip-statistic
    expect_true(length(which(res[,c("intial.p.values"),with=FALSE]>0)) == 2)

    #this method should select "<V450-A>"
    expect_true(res[channel=="<V450-A>"]$area.ratio == 1)

    #it passes over <B710-A> because of negative difference in ICL between 2 and 3 mode
    expect_true(res[channel=="<B710-A>"]$diff.icl < 0)

    #it doesn't compute difference in ICL for <G560-A>
    expect_true(res[channel=="<G560-A>"]$diff.icl== -Inf)

    #it ranks seven channels
    expect_true(min(sort(res[,c("score"),with=FALSE]$score) == seq(7))==1)
}
)
