#' @importFrom ggplot2 autoplot
#' @importFrom ggcyto as.ggplot labs_cyto
best.ICL <- function(flist, fr, debug.mode = FALSE, plotEnv = new.env(parent = emptyenv()), ...){
  res <- sapply(flist, function(gate){
    res <- gate@posteriors[[1]]
    score <- res[["ICL"]]

    gate_pct <- res[["gate_pct"]]
    area.ratio <- ifelse(gate_pct > 0.5, 1 - gate_pct, gate_pct)
    data.table(score, area.ratio)
  }, simplify = FALSE)
  metrics <- rbindlist(res, idcol = "channel")


  if(debug.mode){

    p <- autoplot(fr)
    #save plots and metrics to global hash table
    for(cn in names(p)){
      plotEnv[[cn]] <- as.ggplot(p[[cn]] + labs_cyto("marker"))
      plotEnv[[cn]] <- plotEnv[[cn]] + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    }

  }

  return(metrics)

}
