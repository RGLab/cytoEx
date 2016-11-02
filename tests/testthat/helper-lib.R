require(flowWorkspace)

dataDir <- system.file("extdata",package="flowWorkspaceData")

suppressMessages(gs <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE)))
suppressMessages(gs.singlets <- flowWorkspace::clone(gs, isNew = FALSE, isEmpty = FALSE))
for(node in flowWorkspace::getChildren(gs.singlets, "singlets"))
  flowCore::Rm(node, gs.singlets)

fs <- flowWorkspace::getData(gs.singlets, "singlets")
fr <- fs[[1, use.exprs = FALSE]]
channels <- getfluorescentChannels(fr)
fr.singlet <- fs[[1,channels]] #we now only have fluorescent channels

suppressMessages(gs.cd3 <- flowWorkspace::clone(gs, isNew = FALSE, isEmpty = FALSE))
for(node in flowWorkspace::getChildren(gs.cd3, "singlets/CD3+"))
    flowCore::Rm(node, gs.cd3)
fs.cd3 <- flowWorkspace::getData(gs.cd3, "singlets/CD3+")
fr.cd3 <- fs.cd3[[1,channels]]

suppressMessages(gs.cd8 <- flowWorkspace::clone(gs, isNew = FALSE, isEmpty = FALSE))
for(node in flowWorkspace::getChildren(gs.cd8, "singlets/CD3+/CD8"))
    flowCore::Rm(node, gs.cd8)
fs.cd8 <- flowWorkspace::getData(gs.cd8, "singlets/CD3+/CD8")
fr.cd8 <- fs.cd8[[1,channels]]

#load DC data for testing
suppressMessages(gs.dc.in <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_DC_auto",full = TRUE)))
suppressMessages(gs.dc <- flowWorkspace::clone(gs.dc.in, isNew = FALSE, isEmpty = FALSE))
for (i in seq(length(gs.dc))){
    for(node in flowWorkspace::getChildren(gs.dc[[i]], "Monocytes/Live")){
        flowCore::Rm(node, gs.dc[[i]])
    }
}
fs.dc.mon <- flowWorkspace::getData(gs.dc, "Monocytes/Live")
fr.dc.mon <- fs.dc.mon[[1, use.exprs = FALSE]]
                                        #mon.channels <- cytoEx:::getfluorescentChannels(fr.dc.mon)
mon.channels <- getfluorescentChannels(fr.dc.mon)
fr.dc.mon <- fs.dc.mon[[1,mon.channels]] #we now only have fluorescent channels for DC data


#load bcell data for testing
suppressMessages(gs.bcell.in <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_bcell_auto",full = TRUE)))
suppressMessages(gs.bcell <- flowWorkspace::clone(gs.bcell.in, isNew = FALSE, isEmpty = FALSE))
for (i in seq(length(gs.bcell))){
    for(node in flowWorkspace::getChildren(gs.bcell[[i]], "Live")){
        flowCore::Rm(node, gs.bcell[[i]])
    }
}
fs.bcell <- flowWorkspace::getData(gs.bcell, "Live")
fr.bcell <- fs.bcell[[1, use.exprs = FALSE]]
                                        #bcell.channels <- cytoEx:::getfluorescentChannels(fr.bcell)
bcell.channels <- getfluorescentChannels(fr.bcell)
fr.bcell <- fs.bcell[[1,bcell.channels]] #we now only have fluorescent channels

