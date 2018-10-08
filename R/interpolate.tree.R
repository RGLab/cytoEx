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
    ungated.markers <- all.markers[!sapply(all.markers, function(marker)isGated(marker, path))]
    parents <- path#init the parent with the current path
    for(marker in ungated.markers)
    {

        #find the gate paths that has the target gates
        gate.paths <- find.gates(marker, all.nodes)

        if(length(gate.paths)>0){
          #select the nearest neighour among all the candidates
          dist <- sapply(gate.paths, function(gate.path){
            path.similarity(gate.path, path)
          })
          ind <- sort(dist, decreasing = TRUE)
          #pick one pair of gates from the nearest neighours
          gate.paths <- names(ind[1:2])
          #validity check to see if they are indeed a pair
          sameParent <- dirname(gate.paths[1]) == dirname(gate.paths[2])
          sameTerminal <- setequal(paste0(marker,c("+", "-")), basename(gate.paths))
          if(!sameParent||!sameTerminal)
            stop("The source gates are not paired!")

          #add the gates and update parent paths
          parents <- copy.gates(gs, parents, gate.paths)

        }


    }
  }

}

#' copy the gates to the existing parents
#' @return the new sub tree paths
copy.gates <- function(gs, parents, gate.paths)
{
  sapply(parents, function(parent){
    sapply(gate.paths, function(gate.path) {
      gate <- getGate(gs, gate.path)
      message("copy gate from ", gate.path, " to ", parent)
      add(gs, gate, parent = parent, name = basename(gate.path))
      file.path(parent, basename(gate.path))
    })
  })

}
#' calculate the similarity between two gating paths
#' @param x,y two strings representing the gating paths
#' @examples
#' #the first is smaller is than the second one
#' path.similarity("/LYM/singlets/CD8+/CD45RA+/CCR7+", "/LYM/singlets/CD8+/CD45RA-/CD3+")
#' path.similarity("/LYM/singlets/CD8+/CD45RA+/CCR7+", "/LYM/singlets/CD8+/CD45RA-/CCR7+")
#' path.similarity("/LYM/singlets/CD8+/CD45RA+/CCR7+", "/LYM/singlets/CD8-/CD45RA-/CD4-/HLADR-/CD3+")
path.similarity <- function(x, y)
{
  #convert to two vector
  x <- strsplit(x, split = "/")[[1]][-1]
  y <- strsplit(y, split = "/")[[1]][-1]

  ind <- x%in%y
  #find pos of the first discrepancy node
  pos <- which(!ind)
  if(length(pos)==0)
    stop("identical path in the same tree?", x)
  return(pos[1]-1)
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