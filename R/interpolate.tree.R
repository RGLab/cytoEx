#' add gates to the gating tree so that it is balanced complete tree
#'
#' It is done by copying the gates from the nearest the neighouring branches to the path where the respective
#' marker is missing from the existing path. So that the gating tree is more balanced.
#'
#' @param gs GatingSet
#' @export
interpolate.tree <- function(gs)
{
  all.markers <- markernames(gs)
  all.leaves <- getLeafNode(gs, showHidden = TRUE)
  all.nodes <- getNodes(gs, showHidden = TRUE)
  for(path in all.leaves)
  {
    #check the missing markers
    for(marker in all.markers)
    {
      if(!isGated(marker, path))
      {
        #find the gate paths associated with marker
        gate.paths <- find.gates(marker, all.nodes)

        if(length(gate.paths)>0){
          #select the nearest neighour
          dist <- sapply(gate.paths, function(gate.path){
            path.dist(gate.path, path)
          })
          ind <- which(dist == min(dist))
          if(length(ind)!=2)
            stop("a pair of nearest neighouring paths is expected!")
          gate.paths <- gate.paths[ind]
          #add the gates
          for(gate.path in gate.paths) {
            gate <- getGate(gs, gate.path)
            message("copy gate from ", gate.path, " to ", path)
            add(gs, gate, parent = path, name = basename(gate.path))
          }
        }

      }
    }
  }

}

#' calculate the distance between two gating paths
#' @param x,y two strings representing the gating paths
#' @examples
#' #the first is smaller is than the second one
#' path.dist("/LYM/singlets/CD8+/CD45RA+/CCR7+", "/LYM/singlets/CD8+/CD45RA-/CD3+")
#' path.dist("/LYM/singlets/CD8+/CD45RA+/CCR7+", "/LYM/singlets/CD8-/CD45RA-/CD4-/HLADR-/CD3+")
path.dist <- function(x, y)
{
  #convert to two vector
  x <- strsplit(x, split = "/")[[1]][-1]
  y <- strsplit(y, split = "/")[[1]][-1]

  for()

}

#' find the gate paths that terminates with the given marker
#' @param marker marker name that is associated with the target gate
#' @param all.nodes the candidate gating paths to be searched from
#' @examples
#' \dontrun{
#' find.gates("CD8", getNodes(gs))
#' #expect to return something like "/LYM/singlets/CD8+" and "/LYM/singlets/CD8-"
#' }
#'
find.gates <- function(marker, all.nodes)
{
  all.nodes[grepl(paste0(marker, "[\\+-]$"), all.nodes)]
}