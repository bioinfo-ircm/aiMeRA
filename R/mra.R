#'@title Modular Reponse Analysis function
#'
#'@description Calculation of connectivity coefficients between modules in a biological network
#' using modular response analysis.
#'@usage mra(tab,matp,check=TRUE,Rp=FALSE)
#'@param tab A data.frame containing experimental data in a specific format (see details).
#'@param matp The perturbation matrix. Names of modules (rows) and perturbations (columns) must correspond to names of rows
#'            and columns in tab.
#'@param check Logical. Should the dataset and perutbation matrix be checked for input errors?
#'@param Rp Logical. TRUE if ```tab``` is the calcuated global response matrix.
#'@return A list containing the connectivity map, the local responses matrix, the network responses matrix to perturbations
#'         and the basal line for all modules.
#'@details It assumes that one perturbation must affect only one biological module of the network. This is specified as binary values in
#'the perturbation matrix. It also assumes that row names in data tables are the names of biological modules and column names are the
#'the names of perturbations.
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

mra=function (tab, matp, check = TRUE, Rp=FALSE)
{
  if (check)
    data.setup(tab, matp)
  if (!is.matrix(matp))
    matp = as.matrix(matp)
  lb = colnames(matp)[colSums(matp) == 0]
  if(!Rp)
  {
    glb = tab[rownames(matp), lb]
    gl = tab[rownames(matp), colnames(matp)[colnames(matp)!=lb]]
    Rp = as.matrix(2 * (gl - glb)/(gl + glb))
  }
  else
    Rp=tab[rownames(matp), colnames(matp)[colnames(matp)!=lb]]
  if (ncol(Rp) == 1)
  {
    rownames(Rp) = rownames(matp)
    colnames(Rp) = colnames(matp)[colSums(matp) != 0]
  }
  if (any(rowSums(matp) == 0))
  {
    aux=matrix(0,ncol=sum(rowSums(matp)==0),nrow=nrow(Rp))
    colnames(aux)=rep("bid",sum(rowSums(matp)==0))
    rownames(aux)=rownames(Rp)
    aux[rowSums(matp)==0,]=diag(ncol(aux))
    Rp = cbind(Rp, aux)
    matp = cbind(matp, aux)
  }
  matp = matp[,colSums(matp) != 0]
  rp = solve(diag(diag(matp %*% solve(Rp))))
  rownames(rp) = rownames(Rp)
  colnames(rp) = names(sort(apply(matp, 2, function(x)which(x ==1))))
  r = -rp %*% matp %*% solve(Rp)
  return(list(link_matrix = r, global_matrix = Rp, local_matrix = rp))
}
