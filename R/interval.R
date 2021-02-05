#'@title Confidence interval calculation for MRA connectivity coefficients
#'@description A function that estimates confidence intervals of MRA connectivity coefficients
#'based on replicates. It uses a bootstrap algorithm.
#'@usage interval(tab,mean=0,sd.tab,matp,Rp=FALSE,n=10000)
#'@param tab A data table containing the experimental data in the format output by ```data.setup```.
#'@param mean mean of the normal distribution for creating the noisy matrices.
#'@param sd.tab A data table containing the standard deviation values of replicates for creating the noisy matrices.
#'@param matp A perturbation rule table, with rows corresponding to the MRA modules and columns to perturbations.
#'@param n Number of samples for the bootstrap algorithm.
#'@param Rp Logical. TRUE if ```tab``` is the calcuated global response matrix.
#'@return A list containing the upper and lower values of the confidence interval of each connectivity coefficient.
#'@export
#'@importFrom stats rnorm
#'@examples
#'#Confidence intervals are obtained from the variability within the global
#'#reponse matrices (R). The basal condition is estradiol (E2)
#'#'#We first average the technical replicates of each biological replicate
#'#and then only keep the (averaged) biological replicates for the calculation
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'tec.av=list(data2sdmean(data[1:2])$mean,data2sdmean(data[3:4])$mean,data2sdmean(data[5:6])$mean)
#'data.mean=data2sdmean(tec.av)$mean
#'data.rp=global.matrix(data.mean,"E2")
#'rules=c("E2+siLCoR->LCoR","E2+siRIP140->RIP140","Et->Luciferase","E2->0")
#'matp=read.rules(rules)
#'#The variance of each variable was estimated employing an estimator optimized for a
#'#small sample size from Statistical Process Control theory
#'#(Wheeler and Chambers, 1992; Harter, 1960). The standard deviation was computed for
#'#the global response matrices and stored into the sd.ex table which is included
#'#in the package)
#'interval(data.rp,sd.tab=sd.ex,matp=matp,Rp=TRUE)

interval=function(tab,mean=0,sd.tab,matp,Rp=FALSE,n=10000)
{
  if(!is.matrix(sd.tab))
    sd.tab=as.matrix(sd.tab)
  if(!Rp)
    data.setup(tab,matp)
  genes = rownames(matp)
  lb = colnames(matp)[colSums(matp) == 0]
  pt=colnames(matp)[colnames(matp)!=lb]
  tab=switch(2-Rp,tab[genes,pt],tab[genes,c(pt,lb)])
  sd.tab=switch(2-Rp,sd.tab[genes,pt],sd.tab[genes,c(pt,lb)])
  rlist=list()
  rplist=list()
  for (i in 1:n)
  {
    bruit=matrix(sapply(as.vector(sd.tab),function(x)rnorm(1,mean,x)),ncol=ncol(sd.tab))
    xbruit=tab+bruit
    res = mra(xbruit, matp, check = FALSE, Rp=ifelse(Rp,TRUE,FALSE))
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
