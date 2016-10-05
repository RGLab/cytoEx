getLeafNode <- function(gs, ...){
  nodes <- getNodes(gs, ...)
  isTerminal <- sapply(nodes,function(thisNode){
    length(getChildren(gs,thisNode))==0
  })
  nodes[isTerminal]
}