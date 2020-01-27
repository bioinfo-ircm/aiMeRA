#'@title Confidence interval calculation for MRA connectivity coefficients
#'@description A function that takes in count variability on experimental data replicas and uses the bootstrap algorithm to calculate
#'             a confidence interval for each MRA connectivity coefficient in a biological network.
#'@usage interval(tab,sd.tab,matp,n=10000,nrep=2)
#'@param tab A data frame containing the experimental data in a specific format.
#'@param sd.tab A data frame containing the standard deviation values for each biological module.
#'@param matp A perturbation matrix with rows corresponding to biological modules and columns to perturbations.
#'@param n Number of iterations for the bootstrap algorithm.
#'@param nrep Number of experimental replicas. Default is 2.
#'@return A list containing the upper and lower values of the confidence interval for each connectivity coefficient.
#'@export
#'@importFrom stats quantile
#'@importFrom stats rnorm
#'@examples
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'sd.mean=data2sdmean(data)
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'interval(sd.mean$mean,sd.mean$sd,matp,nrep=6)

interval=function(tab,sd.tab,matp,n=10000,nrep=2)
{
  if(!is.matrix(sd.tab))
    sd.tab=as.matrix(sd.tab)
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
  inter=c(lapply(rlist, function(x) quantile(x,type=5, prob = seq(0, 1, length = 21))[c("5%","95%")]),
          lapply(rplist, function(x) quantile(x,type=5, prob = seq(0, 1, length = 21))[c("5%","95%")]))
  inter=inter[sapply(inter,function(x)x[1]!=x[2])]
  return(inter)
}
