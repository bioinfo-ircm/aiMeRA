#'@title Creates a perturbation matrix from a set of perturbation rules.
#'@description A function that takes a set of perturbation rules as input and creates a perturbation matrix to be used for
#'modular response analysis.
#'@usage read.rules(obj)
#'@param obj Either the path to a perturbation rules file without header or a vector of strings specifying the perturbation rules (see details).
#'@details Perturbations rules are a set of strings specifying the action of perturbations upon modules in a network. In both cases modules
#'and perturbations must have the same names as the names of rows and columns in experimental data tables (names are case sensitive).
#'
#'The rules syntax must be  ```Perturbation->Module``` for a single perturbation,```Perturbation->0``` to specify the basal line and
#'```0->Module``` to specify modules that were not perturbated but for which connectivity coming from other modules can be retrieved.
#'At least two rules and a basal line must be written.
#'
#'
#'@export
#'@importFrom data.table fread
#'@importFrom data.table is.data.table
#'@importFrom data.table setDF
#'@examples
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'read.rules(rules)

read.rules=function(obj)
{
  if(is.character(obj)&length(obj)==1)
    data=data.table::fread(obj,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)$V1
  else if(is.data.frame(obj))
  {
    if(ncol(obj)>1)
      return(message("Only one character rules vector was expected"))
    if(is.data.table(obj))
      obj=setDF(obj)
    data=obj[,1]
  }
  else if (is.vector(obj))
    data=obj
  else
    return(message("Incorrect input format"))
  if(length(data)<2)
    return(message("Not enough perturbations for the MRA analysis"))
  if(!all(grepl("->",data)))
    return(message(paste("incorrect syntax in",data[!grepl("->",data)])))
  if(!any(grepl("*->0$",data)))
    return(message("The basal line is not specified in the perturbation rules file"))
  if(sum(grepl("*->0$",data))>1)
    return(message("The basal line is two times defined "))
  if(any(grepl(" ",data)))
    data=gsub(" ","",data)
  if(any(grepl("^0->*",data)))
    bid=grepl("^0->*",data)
  data=strsplit(data,"->")
  if(any(duplicated(sapply(data,function(x)x[1])[sapply(data,function(x)x[1])!=0])))
    return(message("At least one perturbation affect more than one module in the perturbation rules file"))
  if(any(duplicated(sapply(data,function(x)x[2]))))
    return(message("One or more modules are affected by more than one perturbation in the perturbation rules file"))
  if(exists("bid"))
  {
    mat=matrix(0,nrow=length(data)-(1+sum(bid)),ncol=length(data)-(1+sum(bid)))
    rownames(mat)=sapply(data, function(x)x[2])[sapply(data, function(x)x[2])!="0"&sapply(data, function(x)x[1])!="0"]
    colnames(mat)=sapply(data,function(x)x[1])[sapply(data, function(x)x[2])!="0"&sapply(data, function(x)x[1])!="0"]
    aux=matrix(0,nrow=sum(bid),ncol=ncol(mat))
    rownames(aux)=sapply(data, function(x)x[2])[sapply(data, function(x)x[1])=="0"]
    mat=rbind(mat,aux)
  }
  else
  {
    mat=matrix(0,nrow=length(data)-1,ncol=length(data)-1)
    rownames(mat)=sapply(data, function(x)x[2])[sapply(data, function(x)x[2])!="0"]
    colnames(mat)=sapply(data,function(x)x[1])[sapply(data, function(x)x[2])!="0"]
  }
  diag(mat)=1
  mat=cbind(mat,0)
  colnames(mat)[ncol(mat)]=sapply(data,function(x)x[1])[sapply(data, function(x)x[2])=="0"]
  return(mat)
}
