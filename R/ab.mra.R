#'@title Inference of a double perturbation outcome from a connectivity network created with MRA.
#'@description Function that first calculates the connectivity map of a network containing biological
#'modules and perturbation as specified in the perturbation matrix of the network. Then, it infers the output
#'of combining two perturbations used in the construction of the MRA model.
#'@usage ab.mra(data,matp,pred=NULL,pert1,pert2,inval=c(-1,1),step=0.1,Rp=FALSE,ab=TRUE)
#'@param data A data frame containing the experimental data in the specific format for MRA calculations.
#'@param matp The perturbation matrix. Names of modules (rows) and perturbations (columns) must correspond to names of rows
#'            and columns in tab.
#'@param pred String. Name of the double perturbation to be inferred (if the experimental values for all modules of
#'such perturbation are given as reference in ```data```). Default is NULL.
#'@param pert1 String. Name of the first individual peturbation.
#'@param pert2 String. Name of the second individual perturbation.
#'@param inval A two values vector given the lower and the upper limit of the interval of the a and b coefficients
#' for the ab.mra calculation (See details).
#'@param step Number. Increment of the sequence within the interval for the a and b coefficients.
#'@param Rp Logical. TRUE if ```data``` is the calcuated global response matrix. Default is FALSE
#'@param ab Logical. If TRUE then the inferred values of a double perturbation are obtained
#'by using the a and b coefficients as defined for the ab.mra calculation. If ab=FALSE then
#'a=1 and b=1 (See details).
#'@return List. If ```pert``` is provided and ```ab=TRUE``` then the inferred value, the refrence data
#'and the values for the a and b coefficients are returned. If ```pert``` is provided and ```ab=FALSE``` then the inferred and
#'the reference values are returned. If ```pert``` is not provided then only the inferred values are returned. This
#'will produce a warning message suggesting that the inferred values must be validated by experimental data.
#'@details The ab.MRA inference is based on the hypothesis that it is possible to infer experimental values for
#'the combination of two perturbations used for the classical MRA calculation of the connectivity of a network. The hypothesis
#'is that by including the two local responses of two perturbations into a new system for which the network
#'connectivity is known then it is possible to infer the values of a double perturbation for all the biological
#'modules in the network.
#'
#'Two coefficients (a and b) that ponderate the local responses to perturbations (local_matrix) indicate
#'whether the effect of the combined perturbations remains equal as their individual effect (a=1 and b=1), or it is attenuated (abs(a)<1 and abs(b)<1)
#'or amplified (abs(a)>1 and abs(b)>1). The a, b coefficients are coefficients that minimize the euclidian
#'distance between the infered values and the reference experimental values of all modules.
#'@export
#'@examples
#'#Inferrence of a double perturbation by two siRNAs (siRIP140 and siLCoR)
#'#in a E2 stimulated biological condition of the ERa-RIP140-LCoR network.
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'data.mean=data2sdmean(data)$mean
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'#Inferrence with the experimental value of reference and without the a and b coefficients
#'ab.mra(data.mean,matp=matp,pred="E2+siLCoR+siRIP140",pert1="E2+siLCoR",
#'        pert2="E2+siRIP140",ab=FALSE)
#'#Inferrence with the experimental value of reference and the a and b coefficients
#'ab.mra(data.mean,matp=matp,pred="E2+siLCoR+siRIP140",pert1="E2+siLCoR",
#'        pert2="E2+siRIP140")
#'#Inferrence without the expeimental value of reference and without the a and b coefficients.
#'ab.mra(data.mean,matp=matp,pred="E2+siLCoR+siRIP140",pert1="E2+siLCoR",
#'        pert2="E2+siRIP140")
#'
#'#Inferrence of a biological module (GREB1) without the a,b coefficients,
#'#which does not have an individual perturbation
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'data.mean=data2sdmean(data)$mean
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0","0->GREB1")
#'ab.mra(data.mean,matp=matp,pred="E2+siLCoR+siRIP140",pert1="E2+siLCoR",pert2="E2+siRIP140",
#'       ab=FALSE)



ab.mra=function(data,matp,pred=NULL,pert1,pert2,inval=c(-1,1),step=0.1,Rp=FALSE,ab=TRUE)
{
  if(!pert1%in%colnames(matp)|!pert2%in%colnames(matp))
    stop("pert1 and pert2 should be perturbations as writen in column names of the the perturbation matrix")
  res=mra(data,matp,check=F,Rp=ifelse(Rp,TRUE,FALSE))
  rpa=diag(res$local_matrix)
  coeff=seq(inval[1],inval[2],step)
  lb = colnames(matp)[colSums(matp) == 0]
  if(!is.null(pred))
  {
    if(!pred%in%colnames(data))
      stop("pred should be the name of the column in data to be inferred")
    b=data[rownames(matp),pred]
  }
  if(length(inval)!=2 | !is.numeric(inval))
    stop("The interval must be a vector with only two numerical values")
  mult1=rep(0,ncol(res$local_matrix))
  aux=colnames(res$local_matrix)%in%c(pert1,pert2)
  if(ab&!is.null(pred))
  {
    errtot=vector()
    for(k in coeff)
    {
      err=unlist(lapply(coeff,function(l) {
        mult1[aux]=c(k,l)
        a=solve(res$link_matrix,-rpa*mult1)
        if(!Rp)
        {
          glb=data[rownames(matp),lb]
          a=(-glb*((a+2)/(a-2)))
        }

        err=sqrt(sum((a-b)^2))
      }
      ))
      names(err)=paste0(k,",",coeff)
      errtot=c(errtot,err)
    }
    errtot1=errtot[abs(errtot)==min(abs(errtot))]
    albeta=as.numeric(unlist(strsplit(names(errtot1),",")))
    mult1[aux]=albeta
    pred=solve(res$link_matrix,-rpa*mult1)
    if(!Rp)
    {
      glb=data[rownames(matp),lb]
      pred=(-glb*((pred+2)/(pred-2)))
    }
    all_coeff_total=rep(0,ncol(res$local_matrix))
    all_coeff_total[aux]=albeta
    names(all_coeff_total)=rownames(res$link_matrix)
    return(list(coeff.ab=round(all_coeff_total,digits = 1),inferred=pred,data=b))
  }
  else
  {
    rpa[!aux]=0
    a=solve(res$link_matrix,-rpa)
    if(!Rp)
    {
      glb=data[rownames(matp),lb]
      a=(-glb*((a+2)/(a-2)))
    }
    if(!is.null(pred))
      return(list(inferred=a,data=b))
    else
    {
      warning("The inferred values should be validated by experimental data")
      return(a)
    }
  }
}
