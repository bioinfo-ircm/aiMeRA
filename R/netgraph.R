#'@title Plotting of networks connectivity map
#
#'@description Plots a graphic representation of a biological network connectivity map calculated by MRA.
#'
#'@usage netgraph(map,layout=igraph::layout_with_kk,pertu=NULL,inter=NULL,
#'       cutoff=NULL,main=NULL,digits=2,outfile=NULL,module.in=NULL,
#'       module.out=NULL, no.module=NULL,neg.col="red",pos.col="green", ...)
#'@param map A list containing the connectivity map (link matrix) of a network, and local matrix responses
#'to perturbations, i.e., the output of ```mra```.
#'@param layout Either a data.frame containing the x and y coordinates and the color of the vertices in the network,
#'or a layout function according to the igraph package. Default uses the "layout_with_kk" igraph layout.
#'@param pertu The name of perturbations to be plotted as vertices in the network.
#'@param inter Confidence intervals calculated by ```interval```; connectivity coefficient with a confidence
#'interval that does not include 0 are deemed significant and marked by an asterisk.
#'@param cutoff Minimum value for a connectivity link coefficient between two modules for being plotted.
#'@param main A character string giving the title of the plot.
#'@param digits Number of digits to be plotted for each connectivity coefficient.
#'@param module.in The name of one or more nodes for which only the incoming connectivity
#'              from other nodes of the network should be plotted.
#'@param module.out The name of one or more nodes for which only the outcoming connectivity
#'              to other nodes of the network should be plotted.
#'@param no.module The name of one or more nodes for which incoming and outcoming connectivity
#'              should not be plotted.
#'@param neg.col A color to be used for arrows with negative connectivity coefficients
#'@param pos.col A color to be used for arrows with positive connectivity coefficients
#'@param outfile Optional. A character string giving the name of the output file in graphML format.
#'@param ... Arguments to be passed to igraph plot such as vertex and edges plotting parameters.
#'@return If layout is an igraph function then a data.frame containing the coordinates of the layout selected is returned.
#'If outfile is not NULL, then the network is written in the graphML file format.
#'@export
#'@importFrom igraph V
#'@importFrom igraph graph
#'@importFrom igraph set_edge_attr
#'@importFrom igraph write.graph
#'@importFrom igraph %>%
#'@importFrom graphics plot
#'@examples
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'#We first average the technical replicates of each biological replicate
#'#and then only keep the (averaged) biological replicates for the calculation
#'tec.av=list(data2sdmean(data[1:2])$mean,data2sdmean(data[3:4])$mean,data2sdmean(data[5:6])$mean)
#'sd.mean=data2sdmean(tec.av)
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'map=mra(sd.mean$mean,matp)
#'inter=interval(sd.mean$mean,sd.ex,matp,nrep=6)
#'netgraph(map,inter=inter)

netgraph=function(map,layout=igraph::layout_with_kk,pertu=NULL,inter=NULL,cutoff=NULL,main=NULL,digits=2,
                   outfile=NULL,module.in=NULL,module.out=NULL,no.module=NULL,neg.col="red",pos.col="green", ...)
{
  map=lapply(map,round,digits=digits)
  coeff=numeric()
  noms=expand.grid(list(rownames(map$link_matrix),colnames(map$link_matrix)))
  noms=noms[!apply(noms,1,function(x)any(duplicated(x))),c("Var2","Var1")]
  {
    if(!is.null(module.in)&is.null(module.out))
      noms=noms[noms$Var1%in%module.in,]
    else if (is.null(module.in)&!is.null(module.out))
      noms=noms[noms$Var2%in%module.out,]
    else if (!is.null(module.in)&!is.null(module.out))
      noms=noms[noms$Var2%in%module.out|noms$Var1%in%module.in,]
  }
  if(!is.null(no.module))
    noms=noms[!(noms$Var2%in%no.module|noms$Var1%in%no.module),]
  for(i in 1:nrow(noms))
    coeff=c(coeff,map$link_matrix[noms[i,2],noms[i,1]])
  if(!is.null(cutoff))
  {
    if(cutoff<max(abs(coeff)))
      coeff[abs(coeff)<cutoff]=0.00
    else
      stop("Cutoff too big. At least one network connnexion must exist ")
  }
  if(any(coeff==0.00))
  {
    bid=apply(map$global_matrix,2,function(x)which(x==0.00))
    bid=unlist(lapply(bid,function(x)length(x)))
    noms=noms[coeff!=0.00,]
    coeff=coeff[coeff!=0.00]
  }
  names(coeff)=apply(noms,1,function(x)paste0(x,collapse = "->"))
  noms=as.vector(t(noms))
  if(!is.null(pertu))
  {
    if(!all(pertu%in%colnames(map$local_matrix)))
      stop("Perturbations to plot are not in the map object")
    mols=rownames(map$local_matrix)[colnames(map$local_matrix)%in%pertu]
    loc=diag(map$local_matrix)[colnames(map$local_matrix)%in%pertu]
    names(loc)=paste0(pertu,"->",mols)
    noms=c(noms,rbind(pertu,mols))
    coeff=c(coeff,loc)
  }
  width=log(abs(coeff/min(abs(coeff[coeff!=0]))))
  width=log(abs(coeff/0.01))
  ncols=ifelse(coeff<0,neg.col,pos.col)
  g1 <- igraph::graph(noms,directed = T)
  g1=g1%>%set_edge_attr(name = "coeff",value = coeff )
  if(!is.null(inter))
  {
    bin=lapply(inter[!unlist(lapply(inter,function(x)sum(round(x,digits = 3))==0))],round,digits=2)
    bin=!sapply(bin,function(x)(x[1]<0&x[2]>0)|(x[1]>0&x[2]<0))
    names(bin)=gsub(".5%","",names(bin))
    g1=g1%>%set_edge_attr(name = "sig",value = ifelse(bin,"*","")[1:length(coeff)])
    coeff=round(coeff,digits = digits)
    coeff[names(bin)][bin]=paste0(coeff[names(bin)][bin],"*")
  }
  if(is.function(layout))
  {
    coords=do.call(layout,list(g1))
    plot(g1,edge.label=coeff,layout=coords, ... ,vertex.size=10,edge.width=width,edge.curved=0.2,edge.color=ncols,
         edge.arrow.size=0.5,main=main)
    rownames(coords)=igraph::V(g1)$name
    colnames(coords)=c("X","Y")
    if(!is.null(outfile))
      igraph::write.graph(g1,paste0(outfile,".graphml"),format = "graphml")
    return(invisible(coords))
  }
  plot(g1,edge.label=coeff, ... ,vertex.size=10,edge.width=width,edge.curved=0.2,edge.color=ncols,
       edge.arrow.size=1,layout=layout[igraph::V(g1)$name,c("X","Y")]*.1,rescale=F,edge.label.cex=1,vertex.label.cex=1,
       vertex.color=layout[igraph::V(g1)$name,"Col"],main=main)
  if(!is.null(outfile))
    igraph::write.graph(g1,paste0(outfile,".graphml"),format = "graphml")
}
