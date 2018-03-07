#' Compute Hellinger distance
#'
#' Given a mixture model of multivariate Gaussians compute the Helliger distance between each pair of Gaussians
#'
#' @param model a Mclust object
#'
#' @return a matrix NxN where N is the number of cluster in the model where each entry represents the Hellinger distance between the cluster. Diagonal is equal to 0
#' @noRd
computeHellinger<-function(model)
{
  cluster_probabilities <- model$z;
  num_clusters <- model$G

  cluster_probabilities[which(cluster_probabilities==0)] <- .Machine$double.xmin
  HellingerD=matrix(data=0,nrow=num_clusters,ncol=num_clusters)

  for(i in 1:num_clusters)
  {
    for(j in i:num_clusters)
    {
      if(i==j)
      {
        HellingerD[i,j] <- 0
      }else{
        sigma1 <- model$parameters$variance$sigma[,,i]
        sigma2 <- model$parameters$variance$sigma[,,j]
        mean_sigma <- (sigma1+sigma2)/2
        u <- matrix(model$parameters$mean[,i]-model$parameters$mean[,j], ncol=1)

        HellingerD[i,j] <- sqrt(1-((det(sigma1)^(1/4)*det(sigma2)^(1/4))/sqrt(det(mean_sigma)))*exp(- 1/8 * t(u)%*%solve(mean_sigma,tol=1e-50)%*%u))
        HellingerD[j,i] <- HellingerD[i,j]
      }
    }
  }
  return (HellingerD)
}

#' Compute Battacharyya distance
#'
#' Given a mixture model of multivariate Gaussians compute the Battacharyya distance between each pair of Gaussians. It is an estimation of the overlapping between each pair of Gaussians
#'
#' @param model a Mclust object
#'
#' @return a matrix NxN where N is the number of cluster in the model where each entry represents the Hellinger distance between the cluster. Diagonal is equal to 0
#' @noRd
computeBattacharyya<-function(model)
{
  cluster_probabilities <- model$z;
  num_clusters <- model$G

  cluster_probabilities[which(cluster_probabilities==0)] <- .Machine$double.xmin
  Battacharyya=matrix(data=0,nrow=num_clusters,ncol=num_clusters)

  for(i in 1:num_clusters)
  {
    for(j in i:num_clusters)
    {
      if(i==j)
      {
        Battacharyya[i,j] <- 0
      }else{
        sigma1 <- model$parameters$variance$sigma[,,i]
        sigma2 <- model$parameters$variance$sigma[,,j]
        mean_sigma <- (sigma1+sigma2)/2
        u <- matrix(model$parameters$mean[,i]-model$parameters$mean[,j], ncol=1)
        Battacharyya[i,j] <- (1/8)*t(u)%*%solve(mean_sigma,tol=1e-50)%*%u + 0.5*log(det(mean_sigma)/sqrt(det(sigma1)*det(sigma2)))
        Battacharyya[j,i] <- Battacharyya[i,j]
      }
    }
  }
  return (Battacharyya)
}

