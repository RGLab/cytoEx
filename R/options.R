init.openCyto.exaustive <- function(gsid){

  ops <- getOption("openCyto")
  if(!"exaustive"%in%names(ops))
    ops[["exaustive"]] <- list()
  if(!gsid%in%names(ops[["exaustive"]]))
    ops[["exaustive"]][[gsid]] <- new.env(parent = emptyenv())
  options("openCyto" = ops)

}

set.openCyto.exaustive <- function(gsid, node, winner, plotEnv, metrics){
  if(is.null(node)||is.null(gsid))
    stop("Can't record the marker selection process because gs id and node are not set!")

  #init global option if needed
  init.openCyto.exaustive(gsid)

  #annotate the winner plot
  if(!is.null(winner)) {
    plotEnv[[winner]] <- plotEnv[[winner]] + theme(panel.background = element_rect(fill = "grey50", colour = NA))
  }

  ops <- getOption("openCyto")
  ops[["exaustive"]][[gsid]][[node]] <- list(winner = winner, plotEnv = plotEnv, metrics = metrics)
}