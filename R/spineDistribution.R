dendriteDistribution<-function(model)
{
  idxApi<-which(grepl("api",rownames(model$data)))
  idxBas<-setdiff(1:nrow(model$data),idxApi)
  api<-as.numeric(table(model$classification[idxApi]))
  bas<-as.numeric(table(model$classification[idxBas]))
  return(list(api=api,bas=bas))
}

ageDistribution<-function(model,ageDataset)
{
  ageDataset[,1]<-as.character(ageDataset[,1])
  idxAge<-matrix(0,ncol=1,nrow=nrow(model$data))
  for(i in 1:nrow(ageDataset))
  {
    idxAge[which(grepl(ageDataset[i,1],rownames(model$data))),]<-ageDataset[i,2]
  }

  C40<-as.numeric(table(model$classification[idxAge==1]))
  C85<-as.numeric(table(model$classification[idxAge==2]))

  return(list(C40=C40,C85=C85))
}

bothDistribution<-function(mode,ageDataset)
{
  idxApi<-which(grepl("api",rownames(model$data)))
  idxBas<-setdiff(1:nrow(model$data),idxApi)

  ageDataset[,1]<-as.character(ageDataset[,1])
  idxAge<-matrix(0,ncol=1,nrow=nrow(model$data))
  for(i in 1:nrow(ageDataset))
  {
    idxAge[which(grepl(ageDataset[i,1],rownames(model$data))),]<-ageDataset[i,2]
  }

  apiC40<-as.numeric(table(model$classification[intersect(idxApi,which(idxAge==1))]))
  apiC85<-as.numeric(table(model$classification[intersect(idxApi,which(idxAge==2))]))
  basC40<-as.numeric(table(model$classification[intersect(idxBas,which(idxAge==1))]))
  basC85<-as.numeric(table(model$classification[intersect(idxBas,which(idxAge==2))]))

  return(list(apiC40=apiC40,apiC85=apiC85,basC40=basC40,basC85=basC85))
}
