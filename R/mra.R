#'@title Modular Reponse Analysis function
#'
#'@description Calculation of connectivity coefficients between MRA modules
#' using modular response analysis.
#'@usage mra(tab,matp,check=TRUE)
#'@param tab Data in the format returned by ```data.setup``` or a single data table in the format required by ```data.setup```.
#'@param matp A perturbation rule table, with rows corresponding to the MRA modules and columns to perturbations.
#'@param check Logical. If ```TRUE``` then data and perturbation rule are checked for consistency (by calling ```data.setup```).
#'@return A list containing the connectivity map, the local responses matrix, the network responses matrix to perturbations
#'         and the basal line for all modules.
#'@details A list containing the connectivity map, the local responses matrix, the network responses matrix to perturbations
#'and the basal line for all the modules.
#'@export
#'@examples
#'#It creates the connectivity map between 2 transcriptional nuclear coregulators
#'#(RIP140 and LCoR) and estrogen receptor alpha transcriptional activity reported by
#'#a luciferase gene. q-PCR data is stored in the package files used below in the function data.setup
#'#The model is obtained using the E2 stimulated condition.
#'
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'sd.mean=data2sdmean(data)
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'mra(sd.mean$mean,matp,check=TRUE)

mra=function(tab,matp,check=TRUE)
{
  if(check)
    data.setup(tab,matp)
  if(!is.matrix(matp))
    matp=as.matrix(matp)
  lb=colnames(matp)[colSums(matp)==0]
  gl=tab[rownames(matp),colnames(matp)[colnames(matp)!=lb]]
  glb=tab[rownames(matp),lb]
  Rp=as.matrix(2*(gl-glb)/(gl+glb))
  if(ncol(Rp)==1)
  {
    rownames(Rp)=rownames(matp)
    colnames(Rp)=colnames(matp)[colSums(matp)!=0]
  }
  if(any(rowSums(matp)==0))
  {
    aux=matrix(0,ncol=sum(rowSums(matp)==0),nrow=nrow(Rp))
    colnames(aux)=rep("bid",sum(rowSums(matp)==0))
    rownames(aux)=rownames(Rp)
    aux[rowSums(matp)==0,]=diag(ncol(aux))
    Rp=cbind(Rp,aux)
    matp=cbind(matp,aux)
  }
  matp=matp[,colSums(matp)!=0]
  rp=solve(diag(diag(matp%*%solve(Rp))))
  rownames(rp)=rownames(Rp)
  colnames(rp)=names(sort(apply(matp,2,function(x)which(x==1))))
  r=-rp%*%matp%*%solve(Rp)
  return(list(link_matrix=r,global_matrix=Rp,local_matrix=rp,glb=glb))
}
