#'@title Confidence interval calculation for MRA connectivity coefficients
#'@description A function that estimates confidence intervals of MRA connectivity coefficients
#'based on replicates. It uses a bootstrap algorithm.
#'@usage interval(tab,sd.tab,matp,n=10000,nrep=2)
#'@param tab A data table containing the experimental data in the format output by ```data.setup```.
#'@param sd.tab A data table containing the standard deviation values.
#'@param matp A perturbation rule table, with rows corresponding to the MRA modules and columns to perturbations.
#'@param n Number of samples for the bootstrap algorithm.
#'@param nrep Number of experimental replicas. Default is 2.
#'@return A list containing the upper and lower values of the confidence interval of each connectivity coefficient.
#'@export
#'@importFrom stats quantile
#'@importFrom stats rnorm
#'@examples
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'#We first average the technical replicates of each biological replicate
#'#and then only keep the (averaged) biological replicates for the calculation
#'tec.av=list(data2sdmean(data[1:2])$mean,data2sdmean(data[3:4])$mean,data2sdmean(data[5:6])$mean)
#'sd.mean=data2sdmean(tec.av)
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'#The variance of each variable was estimated employing an estimator optimized for a
#'#small sample size from Statistical Process Control theory (Wheeler and Chambers, 1992; Harter, 1960).
#'interval(sd.mean$mean,sd.ex,matp,nrep=6)

interval=function(tab,sd.tab,matp,n=10000,nrep=2)
{
  if(!is.matrix(sd.tab))
    sd.tab=as.matrix(sd.tab)
  data.setup(tab,matp)
  lb=colnames(matp)[colSums(matp)==0]
  genes=rownames(matp)
  rlist=list()
  rplist=list()
  for (i in 1:n)
  {
    pertus=colnames(matp)
    xbruit=tab+rnorm(ncol(tab)*nrow(tab),0,(sd.tab/sqrt(nrep)))
    res=mra(xbruit,matp,check=FALSE)
    diag(res$link_matrix)=NA
    r=as.vector(res$link_matrix)
    rp=diag(res$local_matrix)
    r=r[!is.na(r)]
    if(i==1)
    {
      rlist=as.list(r)
      rplist=as.list(rp)
      next()
    }
    else
    {
      rlist=lapply(1:length(r),function(j) c(rlist[[j]],r[j]))
      rplist=lapply(1:length(rp),function(j) c(rplist[[j]],rp[j]))
    }
  }
  noms=expand.grid(rep(list(genes), 2))
  noms=noms[!apply(noms,1,function(x)any(duplicated(x))),c("Var2","Var1")]
  names(rlist)=apply(noms,1,function(x)paste0(x,collapse = "->"))
  names(rplist)=paste0(colnames(res$local_matrix),"->",rownames(res$local_matrix))
  inter=c(lapply(rlist, function(x) sort(x)[c(round(n*0.025)+1,round(n*0.975)-1)]),
          lapply(rplist, function(x) sort(x)[c(round(n*0.025)+1,round(n*0.975)-1)]))
  inter=inter[sapply(inter,function(x)x[1]!=x[2])]
  inter=lapply(inter,function(x)round(x,digits = 2))
  return(inter)
}
