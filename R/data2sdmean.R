#'@title Standard deviation and average calculation of a data set.
#'@description This function computes the standard deviation and arithmetic mean values of a data set.
#'@usage data2sdmean(x)
#'@param x A list containing 2 or more data.frames of values.
#'@return A list containing the standard deviation and mean values for each variable in x.
#'@export
#'@examples
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'data2sdmean(data)


data2sdmean=function(x)
{
  if(!is.list(x))
    stop("Input is not a list of data tables")
  m=Reduce("+",x)/length(x)
  sd=sqrt(Reduce("+",lapply(x,function(y)(y-m)^2))/(length(x)-1))
  return(list(sd=sd,mean=m))
}
