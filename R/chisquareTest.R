#' Chisquare test
#'
#' Chisquare test to test independence between cluster distribution and dendritic compartment, age or combination of both
#'
#' @param model a Mclust object
#' @param type a string denoting dendrite compartment ("dendrite"), age ("age") or combination of both ("both")
#' @param ageDatasetCSV a path to a csv file where first column has the name of a dendrite and second column has the binary representation of the dendrite
#'
#' @return p-value of the test
#'
#' @examples
#' chiSquareTest(model,type="age",system.file("extdata", "ageDataset.csv", package = "spineSimulation"))
#'
#' @export
chiSquareTest<-function(model,type=c("dendrite","age","both"), ageDatasetCSV=NULL)
{
  if(class(model) != "Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }
  if(!is.character(type))
  {
    stop("Type must be a character")
  }

  type=tolower(type)

  if(!(type %in% c("dendrite","age","both")))
  {
    stop("Parameter type must be one of the following strings: 'dendrite', 'age' or 'both'")
  }
  if(is.null(ageDatasetCSV) & (type %in% c("age","both")))
  {
    stop("Cannot compute test without the age of the dendrites")
   }

  if(type=="dendrite")
  {
    dendrite<-dendriteDistribution(model)
    chitest<-chisq.test(cbind(dendrite$api,dendrite$bas))
    return(chitest$p.value)
  }
  if(type=="age")
  {
    ageDataset<-read.csv(ageDatasetCSV)
    age<-ageDistribution(model,ageDataset)
    chitest<-chisq.test(cbind(age$C40,age$C85))
    return(chitest$p.value)
  }
  if(type=="both")
  {
    ageDataset<-read.csv(ageDatasetCSV)
    both<-bothDistribution(model,ageDataset)
    chitest<-chisq.test(cbind(both$apiC40,both$apiC85,both$basC40,both$basC85))
    return(chitest$p.value)
  }
}

#' Chisquare test by cluster
#'
#' Chisquare test to test independence between cluster distribution and dendritic compartment, age or combination of both for each cluster
#'
#' @param model a Mclust object
#' @param type a string denoting dendrite compartment ("dendrite"), age ("age") or combination of both ("both")
#' @param ageDatasetCSV a path to a csv file where first column has the name of a dendrite and second column has the binary representation of the dendrite
#'
#' @return p-value of the test
#'
#' @export
chiSquareTestCluster<-function(model,type=c("dendrite","age","both"), ageDatasetCSV=NULL)
{

  if(class(model) != "Mclust")
  {
    stop("Parameter model must be a Mclust object")
  }
  if(!is.character(type))
  {
    stop("Type must be a character")
  }

  type=tolower(type)

  if(!(type %in% c("dendrite","age","both")))
  {
    stop("Parameter type must be one of the following strings: 'dendrite', 'age' or 'both'")
  }
  if(is.null(ageDatasetCSV) & (type %in% c("age","both")))
  {
    stop("Cannot compute test without the age of the dendrites")
  }

  if(type=="dendrite")
  {
    dendrite<-dendriteDistribution(model)
    return(chiSquareClusters(cbind(dendrite$api,dendrite$bas)))
  }
  if(type=="age")
  {
    ageDataset<-read.csv(ageDatasetCSV)
    age<-ageDistribution(model,ageDataset)
    return(chiSquareClusters(cbind(age$C40,age$C85)))
  }
  if(type=="both")
  {
    ageDataset<-read.csv(ageDatasetCSV)
    both<-bothDistribution(model,ageDataset)
    return(chiSquareClusters(cbind(both$apiC40,both$apiC85,both$basC40,both$basC85)))
  }
}

#' Chi Square by cluster
#'
#' Compute Chi Square test for each cluster
#'
#' @param a matrix where each row is a cluster and each column contains the information about the distribution to compare
#'
#' @return a vector where each value is the p-value returned by the test
#' @noRd
chiSquareClusters<-function(matrix=distributionMatrix)
{
  percentage<-colSums(matrix)/sum(matrix)

  numSpineCluster<-rowSums(matrix)
  expectedMatrix<-sapply(percentage,function(x){return(numSpineCluster*x)})
  chiSquareMatrix<-expectedMatrix
  for (i in 1:ncol(expectedMatrix))
  {
    chiSquareMatrix[,i]<-((matrix[,i]-expectedMatrix[,i])^2) / expectedMatrix[,i]
  }

  chiSquareTotal<-rowSums(chiSquareMatrix)
  return(pchisq(chiSquareTotal, df=(ncol(matrix)-1), lower.tail = FALSE))
}
