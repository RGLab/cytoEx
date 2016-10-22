#' @importFrom ggplot2 autoplot
#' @importFrom ggcyto as.ggplot labs_cyto
#' @importFrom flowClust flowClust
best.ICL <- function(fr, debug.mode = FALSE, plotEnv = new.env(parent = emptyenv()), parallel_type, mc.cores, ...){
  channels <- colnames(fr)

  f1 <- function(chnl){
    res <- flowClust(fr, varName = chnl, trans = 0, min.count = -1, max.count = -1, randomStart = 0, K = 1:2)
    icl <- sapply(res, slot, "ICL")
    icl.diff <- diff(icl)

    score <- icl.diff[1]

    if(debug.mode){
      p <- ggflowClust.hist(res[[2]], fr)
    }else
      p <- NULL


    area.ratio <- min(res[[2]]@w)
    list(tbl = data.table(score, area.ratio), plot = p)
  }
  if(parallel_type == "none")
    res <- sapply(channels, f1, simplify = FALSE)
  else
    res <- parallel::mcmapply(FUN = f1, channels, mc.cores = mc.cores, mc.preschedule = FALSE,SIMPLIFY = FALSE, USE.NAMES = TRUE)

  metrics <- sapply(res, function(i)i[["tbl"]], simplify = FALSE)
  metrics <- rbindlist(metrics, idcol = "channel")


  if(debug.mode){

    plotObjs <- sapply(res, function(i)i[["plot"]], simplify = FALSE)
    #save plots and metrics to global hash table
    for(cn in names(plotObjs)){
      p <- plotObjs[[cn]] + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
      p <- p + scale_colour_discrete(guide = FALSE) + ylab("")
      marker <- getChannelMarker(fr, cn)[["desc"]]
      p <- p + xlab(marker)
      plotEnv[[cn]] <- p
    }

  }

  return(metrics)

}
