#' Clustering dendritic spines
#'
#' Read a csv file with the features of the morphology of the dendritic spines. Cluster dendritic spines according to a mixture model of multivariate Gaussians.
#'
#' @param csvSpines string value which indicates the path to the csv file where are saved the morphological features of the spines
#' @param G a vector of integers with the number of clusters to try
#' @param scale boolean value that indicates if the data should be standardized
#'
#' @return Mclust object that contains the mixture of multivariate Gaussians
#'
#' @examples
#' model<-spineClustering(csvSpines = system.file("extdata", "data.csv", package = "spineSimulation"), numClusters = c(1:15), scale = T)
#'
#' @export
spineClustering <- function(csvSpines, numClusters=c(1:15), scale=T)
{
  #Read data
  data <- read.csv(csvSpines)
  clusteringData <- data[,-1] #Remove first column because it contains spine names
  clusteringData<-correctPCA(clusteringData)

  if(scale){
    clusteringData<-scale(clusteringData)
  }

  #Clustering with Mclust
  model <- Mclust(clusteringData,G=numClusters)
  rownames(model$data)<-data[,1]

  model$data<-clusteringData
  model$scaled<-scale

  return(model)
}

#' Plot BIC values
#'
#' For a given model plots the score obtained according to BIC criteria for each number of clusters
#' @param model Mclust object that contains the mixture of multivariate Gaussians
#'
#' @return None
#'
#' @examples
#' plotBIC(model)
#'
#' @export
plotBIC<-function(model)
{
  plot(model,what="BIC")
}

#' Function to rotate ellipses in the case that any of them is pointing to the opposite direction it should be.
#'
#' Function to rotate ellipses in the case that any of them is pointing to the opposite direction it should be.
#'
#' @param data is a data.frame with the data of the spines without the spine names
#'
#' @return a
#'
#' @noRd
correctPCA<-function(data)
{
  data_clustering<-data
  PCA_theta_idx<-grep("inst_Theta_*",colnames(data_clustering))
  PCA_phi_idx<-grep("inst_Phi_*",colnames(data_clustering))

  ###Girar los PCA_theta
  for(i in 1:length(PCA_theta_idx))
  {
    index<-which(data[,PCA_theta_idx[i]]<0)
    data_clustering[,PCA_theta_idx[i]]<-abs(data_clustering[,PCA_theta_idx[i]])
    data_clustering[index,PCA_phi_idx[i]]<-(data_clustering[index,PCA_phi_idx[i]]+pi)%%(2*pi)
  }

  return(data_clustering)
}
