#' Overlapping between pairs of clusters
#'
#' Measure the overlapping between each pair of multivariate Gaussians representing the clusters
#'
#' @param model a Mclust object
#'
#' @return Table representing the overlapping between each pair of clusters
#'
#' @examples
#' computeOverlapping(model)
#'
#' @export
computeOverlapping<-function(model)
{
  overlappingTable<-overlap(model$parameters$pro, t(model$parameters$mean),model$parameters$variance$sigma,eps=1e-09)$OmegaMap
  colnames(overlappingTable)<-paste("Cluster", 1:model$G)
  rownames(overlappingTable)<-paste("Cluster", 1:model$G)
  return(overlappingTable)
}
