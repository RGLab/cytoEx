#clean the marker name(sometime it is in the form of 'antigen isotypecontrol',e.g 'CD38 APC')
cleanMarker <- function(marker){
  strsplit(marker, " ")[[1]][[1]]
}

#' check if the marker has already been gated to avoid gating on the same marker repeately on the same path
#' @param marker marker name to check
#' @param gate.path the gating path that contains the seriers of markers
#' @examples
#'  isGated("CD3 APC", "/boundary/nonDebris/lymph/CD19-/CD3+")
#'  isGated("CD38", "/boundary/nonDebris/lymph/CD19-/CD3+")
isGated <- function(marker, gate.path){

  marker <- cleanMarker(marker)
  gated.markers <- strsplit(gate.path, split = "/")[[1]]
  gated.markers <- gated.markers[-1] #rm the first empty string
  gated.markers <- sub("[\\+\\-]$", "", gated.markers)
  marker %in% gated.markers
}