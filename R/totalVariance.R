#' Barplot total variance (log)
#'
#' Barplot representing the log of the total variance of each cluster
#'
#' @param model a Mclust object
#'
#' @return None
#'
#' @examples
#' plotTotalVariance(model)
#'
#' @export
plotTotalVariance<-function(model)
{
  clusterVariance<-matrix(0,ncol=1,nrow=model$G)
  for(i in 1:model$G)
  {
    clusterVariance[i]<-log(det(model$parameters$variance$sigma[,,i]))
  }

  clusterVariance<-data.frame(var=clusterVariance,cluster=as.factor(1:model$G),random=rep(1,model$G))

  ggplot(clusterVariance, aes(x = random, y = var,fill=cluster)) + geom_bar(stat='identity') + coord_flip() + scale_y_reverse()
}
