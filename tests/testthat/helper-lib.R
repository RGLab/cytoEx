#We have to do flowWorkspace:: explictily because there is a bug introduced in devtools::load_all lately
#which automatically loads helper source for testthat. However helper files are evaled in the package private environment which does not
#have all the namespace imports needed and there is currently no straightforward way to insert these imports to the local environment via library call
dataDir <- system.file("extdata",package="flowWorkspaceData")
suppressMessages(gs.tc <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE)))
fs <- flowWorkspace::getData(gs.tc, "singlets")
fr <- fs[[1, use.exprs = FALSE]]
channels <- getfluorescentChannels(fr)

fr.singlet <- fs[[1,channels]] #we now only have fluorescent channels

fr.cd3 <- flowWorkspace::getData(gs.tc, "singlets/CD3+")[[1,channels]]

fr.cd8 <- flowWorkspace::getData(gs.tc, "singlets/CD3+/CD8")[[1,channels]]

#load DC data for testing
suppressMessages(gs.dc <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_DC_auto",full = TRUE)))
fs.dc.mon <- flowWorkspace::getData(gs.dc, "Monocytes/Live")
fr.dc.mon <- fs.dc.mon[[1, use.exprs = FALSE]]
mon.channels <- getfluorescentChannels(fr.dc.mon)
fr.dc.mon <- fs.dc.mon[[1,mon.channels]] #we now only have fluorescent channels for DC data


#load bcell data for testing
suppressMessages(gs.bcell <- flowWorkspace::load_gs(list.files(dataDir, pattern = "gs_bcell_auto",full = TRUE)))
fs.bcell <- flowWorkspace::getData(gs.bcell, "Live")
fr.bcell <- fs.bcell[[1, use.exprs = FALSE]]
bcell.channels <- getfluorescentChannels(fr.bcell)
fr.bcell <- fs.bcell[[1,bcell.channels]] #we now only have fluorescent channels

#set common options
parallel_type <- "multicore"
mc.cores <- ifelse(interactive(), parallel::detectCores(), 1)
P.ITERS <- 1e2

