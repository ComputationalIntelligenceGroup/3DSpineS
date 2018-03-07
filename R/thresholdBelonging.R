#' Get the numbers of spines with a belonging probability under a threshold
#'
#' Given a threshold probability count the number of spines under that bound
#'
#' @param a Mclust object
#'
#' @return num_spines_over_threshold a table with a probability of belonging lower than a threshold
#'
#' @examples
#' thresholdBelonging(model)
#'
#' @export
thresholdBelonging<-function(model)
{
  class_table<-table(model$classification)
  num_spines_over_threshold<-matrix(0, nrow=6, ncol=model$G)
  for(j in 1:model$G){
    num_spines_over_threshold[1,j] <- class_table[j]-length(which(model$z[,j]>0.99))
    num_spines_over_threshold[2,j] <- class_table[j]-length(which(model$z[,j]>0.9))
    num_spines_over_threshold[3,j] <- class_table[j]-length(which(model$z[,j]>0.8))
    num_spines_over_threshold[4,j] <- class_table[j]-length(which(model$z[,j]>0.7))
    num_spines_over_threshold[5,j] <- class_table[j]-length(which(model$z[,j]>0.6))
    num_spines_over_threshold[6,j] <- class_table[j]-length(which(model$z[,j]>min(apply(model$z,1,max))))
  }

  num_spines_over_threshold <- data.frame(num_spines_over_threshold)
  colnames(num_spines_over_threshold) <- sapply(seq(1,model$G,1),function(x){paste("Cluster",x)})
  rownames(num_spines_over_threshold) <- c(0.99,0.9,0.8,0.7,0.6,format(min(apply(model$z,1,max)),digits=3,nsmall=3))

  return(num_spines_over_threshold)
}
