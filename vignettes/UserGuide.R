## ----eval=F--------------------------------------------------------------
#  install.packages("path/to/file/simulateSpines_0.1.tar.gz", repos = NULL, type="source")

## ----message=FALSE-------------------------------------------------------
library(spineSimulation)

## ----eval=FALSE----------------------------------------------------------
#  spineClustering(csvSpines = system.file("extdata", "spineDataset.csv", package = "spineSimulation"), numClusters = c(2:10), scale = T)

## ----fig.width=7.2,fig.height=7.2----------------------------------------
# Plot the BIC score by the number of cluster and their model name.
# Get the BIC score for each number of clusters
df<-data.frame(clusters=as.numeric(rownames(model$BIC)),BIC=apply(model$BIC,1,function(x){return(max(x,na.rm=T))}))

#Plot the BIC score remove the case where there is only 1 cluster as in the paper
ggplot(data=df[which(df$clusters!=1),], aes(x=clusters, y=BIC, group=1)) +geom_line(size=1)+geom_point(size=2)+labs(x="# of clusters",y="BIC score")

#Show the number of clusters that maximize BIC
print(paste("The number of clusters is",model$G))

## ----fig.width=7.2,fig.height=7.2----------------------------------------
MDS<-computeMDS(model,2)
plotMDS(model,MDS)

## ----eval=F--------------------------------------------------------------
#  computeOverlapping(model)

## ----kable, echo=F-------------------------------------------------------
library(knitr)
kable(computeOverlapping(model))

## ----eval=F--------------------------------------------------------------
#  generateArff('/path/to/folder/',model)

## ----fig.height=5,fig.width=8--------------------------------------------
plotGlobalDistribution(model)
plotDendriticCompartment(model)
plotAge(model,system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
plotCombination(model,system.file("extdata", "ageDataset.csv", package = "spineSimulation"))

## ------------------------------------------------------------------------
chiSquareTest(model,"dendrite")
chiSquareTest(model,"age",system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
chiSquareTest(model,"both",system.file("extdata", "ageDataset.csv", package = "spineSimulation"))

## ------------------------------------------------------------------------
chiSquareTestCluster(model,"dendrite")
chiSquareTestCluster(model,"age",system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
chiSquareTestCluster(model,"both",system.file("extdata", "ageDataset.csv", package = "spineSimulation"))

## ----fig.height=5,fig.width=8, warning=FALSE-----------------------------
plotDistance2Soma(model,system.file("extdata", "distanceSpines2soma.csv", package = "spineSimulation"),6)
chiSquareTestDistance(model,system.file("extdata", "distanceSpines2soma.csv", package = "spineSimulation"),6)

## ----eval=F--------------------------------------------------------------
#  # Simulate 5 spines from the cluster 1 and render the second of them
#  newSpines <- spineSampling(model,nSpines=5,cluster=1,seed=1)
#  mesh <- simulation3Dmesh(newSpines,idx=2,iterations=4)
#  shade3d(mesh,col="red")

