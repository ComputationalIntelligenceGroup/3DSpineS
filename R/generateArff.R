#' Write to Arff file
#'
#' Write data of each cluster in a Arff file to read from Weka.
#'
#' @param dirArff is the path to a folder where Arff files are saved
#' @param model is a Mclust object
#'
#' @return None
#'
#' @export
generateArff<-function(dirArff,model)
{
  if(length(attributes(model$data))==4)
  {
    suma<-attributes(model$data)[[3]]
    cociente<-attributes(model$data)[[4]]

    angular_var<-c(8:13,35:46)
    data<-model$data[,-angular_var]
  }
  dataset<-data.frame(data,cls=model$classification)

  for(i in 1:model$G)
  {
    tempDataset<-dataset
    temptempDataset<-dataset
    tempDataset[which(dataset$cls==i),ncol(dataset)]<-1
    tempDataset[which(dataset$cls!=i),ncol(dataset)]<-2
    temptempDataset<-tempDataset

    temptempDataset$cls<-as.factor(temptempDataset$cls)
    write.arff(temptempDataset,file.path(dirArff,paste0("Cluster ",i,".arff")))
  }
}
