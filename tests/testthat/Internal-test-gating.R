library(data.table)
library(flowWorkspace)
path <- "~/rglab/workspace/analysis/Lyoplate_new/output/gated_data/auto"
gslist <- load_gslist(file.path(path, "gslist-tcell"))
#subselect samples across centers
gslist <- subset(gslist, Replicate == 1 & Sample == 12828)
# gs.tc2 <- rbind2(gslist) # convert from gslist to gs

gs <- gslist@data[[6]]
for(n in getChildren(gs, "lymph"))
  Rm(n, gs)
gs <- flowIncubator::swapChannelMarker(gs)
getNodes(gs)
parallel::mclapply
# lapply(sampleNames(gslist)[1:2], function(sn){
  # lapply(gslist[1:2], function(gs){
      openCyto::gating("lymph", gs
                       , min.count = 2000
                       , min.percent = 0.2
                       # , debug.mode = T
                       # , parallel_type = "none"
                       , parallel_type = parallel_type
                       , mc.cores=mc.cores
                       , gating.function = openCyto::mindensity
                       , marker.selection.function = best.dip.ICL
                       , marker.selection.args = list(P.ITERS=P.ITERS)
                      )
      # }
      # , level = 1
      # , mc.cores = 7
      # )