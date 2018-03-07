#' Compute multidimensional scaling
#'
#' Compute distance between clusters and place instances according to probability to belong to each cluster
#'
#' @param model is a Mclust object
#' @param ndim number of dimension t
#'
#' @return a list where first element is a matrix that represents the centroids of each cluster in a ndim-dimensional space and the second element is a matrix where each instance is placed in a ndim-dimensional space
#'
#' @examples
#' computeMDS(model,2)
#'
#' @export
computeMDS<-function(model, ndim)
{
  hellingerDistance <- computeBattacharyya(model)
  MDS <- cmdscale(hellingerDistance, eig=TRUE, k=ndim)
  coordinates <- model$z %*% MDS$points
  return (list(centers=MDS$points, coords=coordinates))
}

#' Plot result of multidimensional scaling
#'
#' Plot the 2D or 3D multidimensional scaling result.
#'
#' @param model an Mclust object
#' @param MDS a list of matrices obtained from computeMDS
#' @param colors a matrix of rgb colors where each row defines the colors of a cluster
#'
#' @return None
#'
#' @examples
#' MDS<-computeMDS(model,2)
#' plotMDS(model,MDS)
#'
#' @export
plotMDS<-function(model, MDS, colors=NULL)
{
  if(ncol(MDS$centers)>3 | ncol(MDS$centers)<2){
    stop("Too high dimension")
  }

  if(is.null(colors)){
    colors <- t(col2rgb(hue_pal()(nrow(MDS$centers))))
  }else{
    if(nrow(colors)!=model$G | ncol(colors)!=3){
      warning("The number of rows or columns in the colors matrix does not match with the number of clusters so default colors are applied")
      colors <- t(col2rgb(hue_pal()(nrow(MDS$centers))))
    }
  }

  coordColors <- model$z %*% colors
  centColors <- colors

  if(ncol(MDS$centers)==3)
  {
    plot3d(MDS$coords, xlab="", ylab="", zlab="", size=12, col=rgb(coordColors[,1],coordColors[,2],coordColors[,3],maxColorValue=255),axes=F)
    bg3d(col="white")
  }

  if(ncol(MDS$centers)==2)
  {
    minimum<-apply(MDS$coords,2,min)
    maximum<-apply(MDS$coords,2,max)
    plot(MDS$coords, xlab="", ylab="", col=rgb(coordColors[,1],coordColors[,2],coordColors[,3],maxColorValue=255),pch=19,xlim=c(minimum[1]-0.1,maximum[1]+0.1),ylim=c(minimum[2]-0.1,maximum[2]+0.1))
  }
}
