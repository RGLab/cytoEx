#' @import ggplot2
ggflowClust.hist <- function(x, data=NULL, subset=1, include=1:(x@K)
                             , histogram=TRUE
                             , labels=TRUE
                             , main=NULL
                             , col=NULL, pch=20, cex=0.6
                             , ...)
{

  pre.obj <- flowClust:::.hist.flowClust(x = x, data = data, subset = subset, include = include, histogram = histogram, labels = labels, ...)


  den <- flowClust:::flowClust.den(x = pre.obj$data2, obj = x, subset = pre.obj$subset, include = include)
  xdim <- as.vector(x@varNames)
  # df <- fortify(data)

  df <- data.frame(x = pre.obj$data2, y = den)

  p <- ggplot(df, aes(x = x, y = ..density..))
  if (histogram){
    p <- p  + geom_histogram(colour = "grey50", fill = "transparent")
  }
  #fitted density
  p <- p + geom_line(aes(y = den))
  p <- p  + xlab(xdim)+ ylim(pre.obj$ylim)
  if (labels) {

    # if (is.null(col)) {
    #   if (length(include)<=4) col <- c("red", "blue", "green", "black")  else col <- 2:(length(include)+1)
    # } else col<-matrix(col, length(include))
    j <- 0
    # for (k in include){
      this.x <- pre.obj$data
      this.y <- pre.obj$ymin - (pre.obj$ylim[2]-pre.obj$ymin)/100*(j<-j+1)
      this.df <- data.frame(x = this.x, label = factor(x@label))
      p <- p + geom_rug(data = this.df, mapping = aes(x = x, y = -2, color = label), position = "jitter", sides = "b", size=0.1)
    # }
  }
  p
}
