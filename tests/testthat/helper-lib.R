require(flowWorkspace)
dataDir <- system.file("extdata",package="flowWorkspaceData")
suppressMessages(gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE)))
suppressMessages(gs.dip <- clone(gs, isNew = FALSE, isEmpty = FALSE))
for(node in getChildren(gs.dip, "singlets"))
  Rm(node, gs.dip)

fs <- getData(gs.dip, "singlets")
fr <- fs[[1, use.exprs = FALSE]]
channels <- getfluorescentChannels(fr)
fr.singlet <- fs[[1,channels]] #we now only have fluorescent channels
