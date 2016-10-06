init.openCyto.exhaustive <- function(gsid){

  ops <- getOption("openCyto")
  if(!"exhaustive"%in%names(ops))
    ops[["exhaustive"]] <- list()
  if(!gsid%in%names(ops[["exhaustive"]]))
    ops[["exhaustive"]][[gsid]] <- new.env(parent = emptyenv())
  options("openCyto" = ops)

}

set.openCyto.exhaustive <- function(gsid, node, winner, plotEnv, metrics){
  if(is.null(node)||is.null(gsid))
    stop("Can't record the marker selection process because gs id and node are not set!")

  #init global option if needed
  init.openCyto.exhaustive(gsid)

  #annotate the winner plot
  if(!is.null(winner)) {
    plotEnv[[winner]] <- plotEnv[[winner]] + theme(panel.background = element_rect(fill = "grey50", colour = NA))
  }

  ops <- getOption("openCyto")
  ops[["exhaustive"]][[gsid]][[node]] <- list(winner = winner, plotEnv = plotEnv, metrics = metrics)
}
