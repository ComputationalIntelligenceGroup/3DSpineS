#' Barplot distance to soma
#'
#' Barplot representing the distribution of the cluster according to their distance to soma divided in intervals
#'
#' @param model a Mclust object
#' @param distanceCSV path to a CSV file where distances to soma are saved
#' @param numIntervals is the number of intervals
#'
#' @return None
#'
#' @examples
#' plotDistance2Soma(model,distanceCSV=system.file("extdata", "distanceSpines2soma.csv", package = "spineSimulation"),6)
#'
#' @export
plotDistance2Soma<-function(model, distanceCSV, numIntervals)
{

  globalDistribution<-as.numeric(table(model$classification))
  globalDistribution<-globalDistribution/sum(globalDistribution)

  clusterName <- list()
  for(i in 1:model$G)
  {
    clusterName[[i]]<-paste("Cluster",i)
  }
  clusters<-unlist(clusterName)

  distanceMatrix<-distance2somaDistribution(model,distanceCSV,numIntervals)
  intervalos<-c(0,colnames(distanceMatrix))

  interval_names<-list()
  for(i in 1:(length(intervalos)-1))
  {
    interval_names[[i]]<-rep(paste(intervalos[i],intervalos[i+1],sep="-"),model$G)
  }

  distanceMatrix<-sweep(distanceMatrix,2,colSums(distanceMatrix),"/")

  distanceDf<-data.frame(distance=unlist(interval_names),values=c(distanceMatrix),clusters=rep(clusters,length(intervalos)-1),total=rep(globalDistribution,length(intervalos)-1))
  distanceDf$values<-as.numeric(as.character(distanceDf$values))
  distanceDf$distance<-factor(as.character(distanceDf[,1]),levels=unique(distanceDf[,1]))

  ggplot(data=distanceDf,aes(distance,values,fill=as.factor(clusters)))+geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(x=distance, ymax=total, ymin=total), position = "dodge") + ylim(0,0.3) + ylab("Percentage") + xlab(expression(paste("Distance from soma in ",mu))) +
    guides(fill=guide_legend(title=NULL)) + #ggtitle("Distance to soma") +
    theme(plot.title=element_text(lineheight=.8,face="bold",size=20),axis.title.y =element_text(face="bold",size=16),axis.title.x =element_text(size=16),
          axis.text.x =element_text(size=14), axis.text.y =element_text(size=12), legend.text =element_text(size=12))
}

#' Chisquare distance to soma
#'
#' Chisquare test checking if the distribution of clusters is independent of the distance to soma
#'
#' @param model a Mclust object
#' @param distanceCSV path to a CSV file where distances to soma are saved
#' @param numIntervals is the number of intervals
#'
#' @return None
#'
#' @examples
#' chiSquareTestDistance(model,distanceCSV=system.file("extdata", "distanceSpines2soma.csv", package = "spineSimulation"),6)
#'
#' @export
chiSquareTestDistance<-function(model,distanceCSV,numIntervals)
{
  distanceMatrix<-distance2somaDistribution(model,distanceCSV,numIntervals)
  return(chisq.test(distanceMatrix)$p.value)
}

#' Distribution of spines according to their distance from soma
#'
#' Compute distribution of spines according to their distance from soma
#'
#' @param model a Mclust object
#' @param distanceCSV path to a CSV file where distances to soma are saved
#' @param numIntervals is the number of intervals
#'
#' @return None
#'
#' @examples
#' chiSquareTestDistance(model,distanceCSV=system.file("extdata", "distanceSpines2soma.csv", package = "spineSimulation"),6)
#'
#' @noRd
distance2somaDistribution<-function(model,distanceCSV="data/distanceSpines2soma.csv",numIntervals)
{
  distance_matrix<-read.csv(distanceCSV)

  if((max(distance_matrix)%%numIntervals)!=0)
  {
    stop(paste("Parameter numIntervals should be a multiple of",max(distance_matrix)))
  }

  #Make data and distance matrix dendrite names match
  dendrite_names<-colnames(distance_matrix)
  dendrite_names<-gsub("__","-",dendrite_names)
  dendrite_names<-gsub("_"," ",dendrite_names)
  spine_name_split<-matrix(unlist(strsplit(as.character(rownames(model$data)),"_Spine")),ncol=2,byrow=T)

  #Apicales hasta 2904
  spine_name_split[,1]<-gsub("api ","",spine_name_split[,1])
  spine_name_split[,1]<-gsub(" api","",spine_name_split[,1])
  spine_name_split[,1]<-gsub("api","",spine_name_split[,1])
  spine_name_split[,2]<-gsub(".mat","",spine_name_split[,2])

  dendrite_names[3]="if6 1 9-1"
  dendrite_names[36]="if6 2 0-2"
  dendrite_names[46]="if6 2 20-3"
  dendrite_names[54]="m16 1 2-1"
  dendrite_names[55]="m16 1 2-2"
  dendrite_names[56]="m16 1 2-3"

  #Apicales hasta 16
  dendrite_names<-gsub("api ","",dendrite_names)

  #Matrix of interval distance
  distancias <- matrix(0,nrow=nrow(model$data),ncol=1)
  for(i in 1:length(spine_name_split[,2])){
    idx <- which(spine_name_split[i,1]==dendrite_names)
    if(length(idx)>1)
    {
      if(i<length(grep("api",rownames(model$data)))){#Las apicales son las que van de 1 a 2904
        idx<-idx[which(idx<17)]#Los nombres de dendrita apical van de 1 a 16
      }else{
        idx<-idx[which(idx>16)]
      }
    }

    distancias[i]<-distance_matrix[as.numeric(spine_name_split[i,2]),idx]
  }

  distancias<-unlist(distancias)
  distancias[1:length(grep("api",rownames(model$data)))]<-distancias[1:length(grep("api",rownames(model$data)))]+100 #Sumo 100 a todas las apicales porque todas empiezan en el tramo 2.

  #Each row provides information about the distance and the cluster
  class_distance<-cbind(distancias,model$classification)

  #Group distances by classification label in a list
  distanceMatrix<-list()
  for (i in 1:model$G)
  {
    distanceMatrix[[i]]<-class_distance[which(class_distance[,2]==i),1]
  }

  #Get minimum value and maximum value
  minimum<-min(unlist(distanceMatrix))-10
  maximum<-max(unlist(distanceMatrix))

  #Generate breaks for the hist plot
  intervalos<-seq(minimum,to=maximum, by=(maximum/numIntervals))

  #For each cluster compute the histogram and return the number of occurence at each distance interval
  distribution<-lapply(distanceMatrix,function(x){return(hist(x,breaks=intervalos,plot=F)$count)})

  #List to matrix conversion
  distributionMatrix<-do.call(rbind, distribution)

  colnames(distributionMatrix)<-intervalos[2:length(intervalos)]
  return(distributionMatrix)
}
