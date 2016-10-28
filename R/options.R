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
  # if(!is.null(winner)) {
    # plotEnv[[winner]] <- plotEnv[[winner]] + theme(panel.background = element_rect(fill = "grey50", colour = NA))
  # }

  ops <- getOption("openCyto")
  ops[["exhaustive"]][[gsid]][[node]] <- list(winner = winner, plotEnv = plotEnv, metrics = metrics)
}

#' get the enviornment that stores the marker selection results
#'
#' @param gs GatingSet
#' @param node the unique path to the gated population node
get.openCyto.exhaustive <- function(gs, node){
  gsEnv <- getOption("openCyto")[["exhaustive"]][[gs@guid]]
  fullpath <- file.path(getParent(gs, node), basename(node))
  gsEnv[[fullpath]]

}

#' plot the flowClust fitting results for a particular gated node
#'
#' @param ... arguments
#'    gs GatingSet
#'    node the unique path to the gated population node
#' @export
plot.flowClust <- function(...){

  marker.selection.res <- get.openCyto.exhaustive(...)
  if(is.null(marker.selection.res))
    return(NULL)
  else{
    plotEnv <- marker.selection.res[["plotEnv"]]
    p <- as.list(plotEnv)
    tb <- marker.selection.res[["metrics"]]
    channel <- tb[order(score, decreasing = TRUE), channel]
    p <- gridExtra::arrangeGrob(grobs = p[channel])
    plot(p)
  }
}

#' get the metrics table computed by the marker selection method
#' @inheritParams plot.flowClust
#' @export
getMetrics <- function(...){
  marker.selection.res <- get.openCyto.exhaustive(...)
  if(is.null(marker.selection.res))
    return(NULL)
  else{
    tb <- marker.selection.res[["metrics"]]
    tb[order(score, decreasing = TRUE),]
  }

}