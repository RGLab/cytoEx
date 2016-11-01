require(flowWorkspace)
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
gs.dip <- clone(gs, isNew = FALSE, isEmpty = FALSE)
for(node in getChildren(gs.dip, "singlets"))
  Rm(node, gs.dip)
fs <- getData(gs.dip, "singlets")
fr <- fs[[1, use.exprs = FALSE]]
pd <- pData(parameters(fr))
pd <- pd[!is.na(pd[["desc"]]),]
channels <- pd[["name"]]
fr.singlet <- fs[[1,channels]] #we now only have fluorescent channels
