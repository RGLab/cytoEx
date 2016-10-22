#' @include utils.R options.R
NULL

#' automated gating
#'
#' The idea is pretty simple: start with the marker with the best separation between -/+ peaks,
#' then iterate over all other markers until you can't find a good separation between peaks.
#' Different methods can be used for determining good/bad separation, such as the ICL or the
#' heuristic method based on peak distance.
#'
#' @param x a GatingSet
#' @param y either missing or a character specifying the starting node
#' @importFrom openCyto gating
#' @export
#' @examples
#' \dontrun{
#' dataDir <- system.file("extdata",package="flowWorkspaceData")
#' gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
#'
#' #clone the gating tree strucuture of the existing manually gated gs(through flowJo)
#' gs2 <- clone(gs, isNew = FALSE, isEmpty = FALSE)
#'
#' #erase all the acestor nodes of CD3+
#' Rm("CD3+", gs2)
#'
#' # run the exhaustive gating from "singlets" node
#' gating("singlets", gs2, min.count = 2000, min.percent = 0.3)
#'
#' # or proceed from any leaf nodes of the existing gating tree
#' gating(gs2, debug.mode = TRUE) #turn on the debug mode to generate plots and messages (particularly in knitr report) for troubleshooting
#'
#' }
#' @rdname gating
#' @import flowWorkspace
setMethod("gating", signature = c("GatingSet","missing"),
          definition = function(x, y, ...) {

            nodes <- getLeafNode(gs)
            for(node in nodes){
              gating(node, gs, ...)
            }
          })


setMethod("gating", signature = c("character", "GatingSet"),
          definition = function(x, y, ...) {
            gating.subnode(x,y,...)
          })

#' @rdname gating
#' @param marker.selection a function that selectes the best marker that has the best bi-module separation based on the various peak statistics (e.g. peak area ratio, peak distance, peak/valley ratio, etc...)
#' @param gating.function the 1d gating function used for exhaustive gating
#' @param min.count the minimum number of cells that allows the gating proceed further
#' @param max.depth the maximum depths of gating path. Default is -1, which is no limits.
#' @param ... other arguments passed to \link{density} function.
gating.subnode <- function(parent, gs
                           , gating.function = mindensity
                           , marker.selection = best.separation
                           , min.count = 1000, min.percent = 0.2
                           , debug.mode = FALSE, max.depth = -1
                           , mc.cores = getOption("mc.cores", 2L)
                           , parallel_type = c("multicore", "cluster", "none"), cl = NULL
                           , ...){
  parallel_type <- match.arg(parallel_type)
  parent <- file.path(getParent(gs, parent), basename(parent))#ensure to get full path
  depths <- length(strsplit(parent, "/")[[1]])
  if(max.depth >0 && depths >= max.depth)
    message("stop gating at ", parent, ". Reaching the maximum gating depths: ", max.depth)
  else
  {
    fs <- getData(gs, parent)
    #for now we use the first sample
    fr <- fs[[1, use.exprs = FALSE]]
    #exclude the non-stained channels
    pd <- pData(parameters(fr))
    pd <- pd[!is.na(pd[["desc"]]),]
    channels <- pd[["name"]]
    fr <- fs[[1,channels]]

    nCell <- nrow(fr)
    if(nCell > min.count){#TODO: use options("openCyto")[[1]][["gating"]][["minEvents"]]

      #check if the markers have already been gated to
      #avoid gating on the same marker repeately on the same path
      is.gated <- sapply(channels, function(channel){
        marker <- getChannelMarker(fr, channel)[, "desc"]
        gated.markers <- strsplit(parent, split = "/")[[1]]
        gated.markers <- gated.markers[-1] #rm the first empty string
        matched <- sapply(gated.markers, function(i){
          i <- sub("[\\+\\-]$", "", i)
          grepl(i, marker)
        })
        any(matched)

      })
      channels <- channels[!is.gated]

      if(length(channels) > 0){
        message("parent: ", parent)
        #get measurements for the cutpoints
        plotEnv <- new.env(parent = emptyenv())
        metrics <- marker.selection(fr, debug.mode = debug.mode, plotEnv = plotEnv
                                    , mc.cores = mc.cores
                                    , parallel_type = parallel_type)
        # add marker column for the purpose of visualization
        metrics[, marker := getChannelMarker(fr, channel)[["desc"]], by = channel]
        metrics.thresholded <- metrics[area.ratio >= min.percent, ]


        if(nrow(metrics.thresholded)==0){
          chnl.selected <- NULL
        }else{
          ind <- which.max(metrics.thresholded[, score])
          chnl.selected <- metrics.thresholded[ind, channel]
          marker <- metrics.thresholded[ind, marker]
        }
        if(debug.mode){
          set.openCyto.exhaustive(gs@guid, parent, winner = chnl.selected, plotEnv = plotEnv, metrics = metrics)
          }

        if(!is.null(chnl.selected)){

          #clean the marker name(sometime it is in the form of 'antigen isotypecontrol',e.g 'CD38 APC')
          marker <- strsplit(marker, " ")[[1]][[1]]
          message("selected marker: ", marker)
          #add the gates and move on the children of the new node
          for(sub in c("+", "-")){
            negated <- sub == "-"
            child <- paste0(marker, sub)
            gate <- gating.function(fr, chnl.selected, ...)
            suppressMessages(add(gs, gate, name = child, parent = parent, negated = negated))
            node.path <- file.path(parent,child)
            suppressMessages(recompute(gs, node.path))
            #resursively gate the sub nodes
            gating(node.path, gs
                   , gating.function = gating.function
                   , marker.selection = marker.selection
                   , min.count = min.count
                   , min.percent = min.percent
                   , debug.mode = debug.mode
                   , max.depth = max.depth
                   , mc.cores = mc.cores
                   , parallel_type = parallel_type
                   , cl = cl
                   , ...)
          }
        }else
          message("skip node '", parent, " due to low pop percent: ", metrics[, max(area.ratio)])
      }

    }else{

        message("skip node '", parent, "' due to the low cell count : ", nCell)
    }
  }

}
