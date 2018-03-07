#' Barplot spine distribution by cluster
#'
#' Barplot of the distributions of the spines by cluster according to a given model
#'
#' @param model a Mclust object
#'
#' @return None
#'
#' @examples
#' plotGlobalDistribution(model)
#'
#' @export
plotGlobalDistribution<-function(model)
{
  if(class(model)!="Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }

  clusterName<-list()

  for(i in 1:model$G)
  {
    clusterName[[i]]<-paste("Cluster",i)
  }

  clusters<-unlist(clusterName)

  distribution<-as.numeric(table(model$classification))
  distribution<-distribution/sum(distribution)

  dfDistribution<-data.frame(values=distribution,clusters=clusters)
  dfDistribution$values<-as.numeric(as.character(distribution))

  ggplot(dfDistribution,aes(clusters,values,fill=as.factor(clusters))) +
    geom_bar(stat="identity")+#ggtitle("Classes of spines") +xlab("Global cluster distribution") +
    ylab("Percentage") + ylim(0,0.3) + guides(fill=guide_legend(title=NULL)) +
    theme(plot.title=element_text(lineheight=.8,face="bold",size=20), axis.title.x =element_text(face="bold",size=16),
          axis.title.y =element_text(face="bold",size=16), axis.text =element_text(size=12),
          legend.text =element_text(size=12))
}


#' Barplot spine distribution by dendritic compartment
#'
#' Barplot of the distributions of the spines by their dendritic compartment and cluster
#'
#' @param model a Mclust object
#'
#' @return None
#'
#' @examples
#' plotDendriticCompartment(model)
#'
#' @export
plotDendriticCompartment<-function(model)
{
  if(class(model)!="Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }

  clusterName<-list()

  for(i in 1:model$G)
  {
    clusterName[[i]]<-paste("Cluster",i)
  }

  clusters<-unlist(clusterName)

  dendrite<-dendriteDistribution(model)
  for(i in 1:length(dendrite))
  {
    dendrite[[i]]<-dendrite[[i]]/sum(dendrite[[i]])
  }

  distribution<-as.numeric(table(model$classification))
  distribution<-distribution/sum(distribution)

  dendriteDf<-data.frame(dendrite=c(rep("Apical",model$G),rep("Basal",model$G)),values=c(dendrite$api,dendrite$bas),clusters=c(clusters,clusters),total=rep(distribution,2))
  dendriteDf$values<-as.numeric(as.character(dendriteDf$values))


  ggplot(data=dendriteDf,aes(dendrite,values,fill=as.factor(clusters)))+geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(x=dendrite, ymax=total, ymin=total), position = "dodge") + ylim(0,0.3) + ylab("Percentage") + xlab("") +
    guides(fill=guide_legend(title=NULL)) +
    theme(plot.title=element_text(lineheight=.8,face="bold",size=20),axis.title.y =element_text(face="bold",size=16),
          axis.text.x =element_text(size=14), axis.text.y =element_text(size=12), legend.text =element_text(size=12))
}

#' Barplot spine distribution by age
#'
#' Barplot of the distributions of the spines by their age and cluster
#'
#' @param model a Mclust object
#' @param ageDatasetCSV a path to a csv file where first column has the name of a dendrite and second column has the binary representation of the dendrite
#'
#' @return None
#'
#' @examples
#' plotAge(model,system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
#'
#' @export
plotAge<-function(model, ageDatasetCSV)
{
  ageDataset<-read.csv(ageDatasetCSV)
  if(class(model) != "Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }

  clusterName <- list()

  for(i in 1:model$G)
  {
    clusterName[[i]]<-paste("Cluster",i)
  }

  clusters<-unlist(clusterName)

  age<-ageDistribution(model,ageDataset)
  for(i in 1:length(age))
  {
    age[[i]]<-age[[i]]/sum(age[[i]])
  }

  distribution<-as.numeric(table(model$classification))
  distribution<-distribution/sum(distribution)

  ageDf<-data.frame(age=c(rep("C40",model$G),rep("C85",model$G)),values=c(age$C40,age$C85),clusters=c(clusters,clusters),total=rep(distribution,2))
  ageDf$values<-as.numeric(as.character(ageDf$values))


  ggplot(data=ageDf,aes(age,values,fill=as.factor(clusters)))+geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(x=age, ymax=total, ymin=total), position = "dodge") + ylim(0,0.3) + ylab("Percentage") + xlab("") +
    guides(fill=guide_legend(title=NULL)) +
    theme(plot.title=element_text(lineheight=.8,face="bold",size=20),axis.title.y =element_text(face="bold",size=16),
          axis.text.x =element_text(size=14), axis.text.y =element_text(size=12), legend.text =element_text(size=12))
}



#' Barplot spine distribution by combination of dendritic compartment and age
#'
#' Barplot of the distributions of the spines by combination of dendritic compartment and age
#'
#' @param model a Mclust object
#' @param ageDatasetCSV a path to the csv where is saved the age of each dendrite
#'
#' @return None
#'
#' @examples
#' plotCombination(model,system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
#'
#' @export
plotCombination<-function(model, ageDatasetCSV)
{
  ageDataset<-read.csv(ageDatasetCSV)
  if(class(model) != "Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }

  ageDataset[,1]<-as.character(ageDataset[,1])
  spineAge<-matrix(0,ncol=1,nrow=nrow(model$data))
  for(i in 1:nrow(ageDataset))
  {
    spineAge[which(grepl(ageDataset[i,1],rownames(model$data))),]<-ageDataset[i,2]
  }

  clusterName <- list()

  for(i in 1:model$G)
  {
    clusterName[[i]]<-paste("Cluster",i)
  }

  clusters<-unlist(clusterName)

  both<-bothDistribution(model,ageDataset)
  for(i in 1:length(both))
  {
    both[[i]]<-both[[i]]/sum(both[[i]])
  }
  distribution<-as.numeric(table(model$classification))
  distribution<-distribution/sum(distribution)

  combinationDf<-data.frame(dendrite=c(rep("apiC40",model$G),rep("apiC85",model$G),rep("basC40",model$G),rep("basC85",model$G)),values=c(both$apiC40,both$apiC85,both$basC40,both$basC85),clusters=c(clusters,clusters,clusters,clusters),total=rep(distribution,4))
  combinationDf$values<-as.numeric(as.character(combinationDf$values))


  ggplot(data=combinationDf,aes(dendrite,values,fill=as.factor(clusters)))+geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(x=dendrite, ymax=total, ymin=total), position = "dodge") + ylim(0,0.3) + ylab("Percentage") + xlab("") +
    guides(fill=guide_legend(title=NULL)) +
    theme(plot.title=element_text(lineheight=.8,face="bold",size=20),axis.title.y =element_text(face="bold",size=16),
          axis.text.x =element_text(size=14), axis.text.y =element_text(size=12), legend.text =element_text(size=12))
}



