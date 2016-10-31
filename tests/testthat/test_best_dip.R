#library(cytoEx)
library(flowWorkspace)
#library(diptest)
#library(parallel)
#library(flowClust)
#library(ggcyto)


test_that("Test best.dip singlet selection unchanged: ", {
    set.seed(172384)
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
    res <- cytoEx:::best.dip(fr,parallel_type="multicore",mc.cores=4)

    #this data set should have two channels with non-zero p-value for dip-statistic
    expect_true(length(which(res[,c("intial.p.values"),with=FALSE]>0)) == 2)
    #data set should select "<B710-A>"
    expect_true(res[channel=="<B710-A>"]$area.ratio == 1)
})



