#'@title Data setup for Modular Response Analysis
#'@description Checks for possible errors in tables of experimental data and in the perturbation matrix.
#'Returns data properly formated for the mra function.
#'@param obj A character string specifying the path to a folder that contains the experimental data tables
#'and the perturbation rules, or a list containing several data tables, or a single data table.
#'@param pert.tab Optional. The perturbation rule table if not provided in ```obj```.
#'@param dec The decimal separator as in ```base::read.csv``` (default is ".")
#'@return If more than one data table are checked, then a list of matrices in the correct
#'format for the mra function is returned. If a single table is checked, then a single matrix
#'in the right format is returned. It also checks for correspondance with the perurbation rules table.
#'@export
#'@importFrom data.table fread
#'@examples
#'#Example using q-PCR data stored in the package files
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'#Examples using a gene network from Kholodenko et.al. PNAS 2002 figure a and b
#'#Here data was obtained from the given R matrix
#'data.A=data.setup(gene.network.A)
#'data.B=data.setup(gene.network.B)

data.setup=function(obj,pert.tab=NULL,dec=".")
{
  if(!is.list(obj)&is.vector(obj))
  {
    if(length(obj)>1)
      return(message("Only one path to the data tables and rules table is allowed"))
    temp=list.files(obj)
    if (any(length(temp)==c(0,1)))
      return(message("At least one data table and a perturbation rules file are required"))
    data=lapply(paste0(obj,temp),data.table::fread,data.table=FALSE,dec=dec)
    if(any(sapply(data,ncol)==1))
    {
      if(sum(sapply(data,ncol)==1)>1)
        return(message("Only one rules file is allowed"))
      if(all(sapply(unlist(data[sapply(data,ncol)==1]),is.character)))
      {
        pert.tab= c(colnames(data[sapply(data,ncol)==1][[1]]),unlist(data[sapply(data,ncol)==1]))
        data=data[sapply(data,ncol)!=1]
      }
      else
        return(message("Wrong format in file with perturbation rules"))
    }
    else
      return(message("File with perturbation rules was not found"))
  }
  else if (is.data.frame(obj)|is.matrix(obj)|is.data.table(obj))
    data=list(obj)
  else
    data=obj
  if (any(sapply(data,ncol)<4))
    return(message("Columns in data tables must contain at least two perturbations for two modules, modules names, and the basal line"))
  if (!any(sapply(data,function(x)is.character(x[,1])|is.character(colnames(x)))))
    return(message("First column of each data table must contain the modules name"))

  if(all(sapply(data,function(x)is.character(x[,1]))))
  {
    data=lapply(data,function(x){row.names(x)=as.character(x[,1]);x})
    data=lapply(data,function(x)x[,-1])
  }
  if(!is.null(pert.tab))
  {
    if((is.data.frame(pert.tab)|is.matrix(pert.tab)))
    {
      temp=sapply(colnames(pert.tab),function(x)paste(x,rownames(pert.tab)[!!pert.tab[,x]],sep = "->"))
      temp[colSums(pert.tab)==0]=paste0(temp[colSums(pert.tab)==0],"0")
      pert.tab=temp
    }
    if(any(grepl(" ",pert.tab)))
      pert.tab=gsub(" ","",pert.tab)
    aux=strsplit(pert.tab,"->")
    if(!any(grepl("*->0",pert.tab)))
      return(message("The basal line is not specified in the perturbation rules file"))
    if(sum(grepl("*->0$",pert.tab))>1)
      return(message("The basal line is two times defined "))
    if(any(duplicated(sapply(aux,function(x)x[1])[sapply(aux,function(x)x[1])!="0"])))
      return(message("At least one perturbation affect more than one module in the perturbation rules file"))
    if(any(duplicated(sapply(aux,function(x)x[2]))))
      return(message("One or more modules are affected by more than one perturbation in the perturbation rules file"))
    if(!all(sapply(aux[!grepl("^0->*",pert.tab)], function(x)x[1]%in%colnames(data[[1]]))))
      return(message("At least one name of perturbations in the perturbations rules file do not correspond to perturbations names in at least one data table"))
    if(!all(sapply(aux[!grepl("*->0$",pert.tab)], function(x)x[2]%in%switch(2-is.character(data[[1]][,1]),data[[1]][,1],rownames(data[[1]])))))
      return(message("At least one name of modules in the perturbation rules file do not correspond to modules names in data tables"))
  }
  if(length(data)>1)
  {
    aux=Reduce(paste,lapply(data, colnames))
    aux=strsplit(aux," ")
    if(!all(sapply(aux, function(x)all(x[1]==x))))
      return(message("At least the name of one column is not the same between data tables"))
    aux=Reduce(paste,lapply(data, rownames))
    aux=strsplit(aux," ")
    if(!all(sapply(aux, function(x)all(x[1]==x))))
      return(message("At least the name of one row is not the same between data tables"))
    data=lapply(data,as.matrix)
    return(data)
  }
  return(as.matrix(data[[1]]))
}
