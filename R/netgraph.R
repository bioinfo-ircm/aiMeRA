#'@title Plotting of networks connectivity map
#
#'@description Plots a graphic representation of a biological network connectivity map calculated by MRA.
#'
#'@usage netgraph(map,layout=igraph::layout_with_kk,pertu=FALSE,inter=NULL,
#'       cutoff=NULL,main=NULL,digits=2,v.type="circle",p.type="rectangle",
#'       outfile=NULL, ...)
#'@param map A list containing the connectivity map (link matrix) of a network, and local matrix responses to perturbations
#'@param layout Either a data.frame containing the x and y coordinates and the color of the vertex in the network,
#'              or a layout function from the igraph package. Default uses the "layout_with_kk" igraph layout.
#'@param pertu A boolean specifying if perturbations should be plotted as vertex in the network.
#'@param inter It takes the confidence interval calculated by the interval function and marks each connectivity coefficient
#'             with an asterisk if the coefficient is significant.
#'@param cutoff Minimum value for a connectivity link coefficient between two modules for being plotted.
#'@param main A character string giving the title of the plot.
#'@param digits Number of digits to be plotted for each connectivity coefficient.
#'@param v.type Vertex type for biological modules. One of “none”, “circle”, “square”, “csquare”, “rectangle”
#'              “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”.
#'@param p.type Vertex type for perturbations. One of “none”, “circle”, “square”, “csquare”, “rectangle”
#'              “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”.
#'@param outfile Optional. A character string giving the name of the output file in xml format.
#'@param ... Arguments to be passed to igraph plot such as vertex and edges plotting parameters.
#'@return If layout is an igraph function it returns a data.frame containing the coordinates of the layout selected.
#'If outfile is not NULL exports the network graph in xml file format.
#'@export
#'@importFrom igraph V
#'@importFrom igraph graph
#'@importFrom igraph set_edge_attr
#'@importFrom igraph write.graph
#'@importFrom igraph %>%
#'@importFrom graphics plot
#'@examples
#'data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
#'sd.mean=data2sdmean(data)
#'rules=c("Et->Luciferase","E2+siRIP140->RIP140","E2+siLCoR->LCoR","E2->0")
#'matp=read.rules(rules)
#'map=mra(sd.mean$mean,matp)
#'inter=interval(sd.mean$mean,sd.mean$sd,matp,nrep=6)
#'netgraph(map,inter=inter)


netgraph=function(map,layout=igraph::layout_with_kk,pertu=FALSE,inter=NULL,cutoff=NULL,main=NULL,digits=2,v.type="circle",
                  p.type="rectangle",outfile=NULL, ...)
{
  map=lapply(map,round,digits=digits)
  coeff=numeric()
  noms=expand.grid(rep(list(rownames(map$link_matrix)), 2))
  noms=noms[!apply(noms,1,function(x)any(duplicated(x))),c("Var2","Var1")]
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
  v.type=rep(v.type,nrow(map$link_matrix))
  if(pertu)
  {
    pertu=colnames(map$local_matrix)
    mols=rownames(map$local_matrix)
    loc=diag(map$local_matrix)
    names(loc)=paste0(pertu,"->",mols)
    if(any(pertu=="bid"))
    {
      loc=loc[pertu!="bid"]
      mols=mols[pertu!="bid"]
      pertu=pertu[pertu!="bid"]
    }
    noms=c(noms,rbind(pertu,mols))
    coeff=c(coeff,loc)
    v.type=c(v.type,rep(p.type,length(pertu)))
  }
  width=abs(coeff/max(abs(coeff)))*5
  ncols=ifelse(coeff<0,"red","green")
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
         vertex.shape=v.type,edge.arrow.size=0.5)
    rownames(coords)=igraph::V(g1)$name
    colnames(coords)=c("X","Y")
    if(!is.null(outfile))
      igraph::write.graph(g1,paste0(outfile,".xml"),format = "graphml")
    return(invisible(coords))
  }
  plot(g1,edge.label=coeff,vertex.shape=v.type, ... ,vertex.size=10,edge.width=width,edge.curved=0.2,edge.color=ncols,
       edge.arrow.size=1,layout=layout[igraph::V(g1)$name,c("X","Y")]*.1,rescale=F,edge.label.cex=1,vertex.label.cex=1,
       vertex.color=layout[igraph::V(g1)$name,"Col"],main=main)
  if(!is.null(outfile))
    igraph::write.graph(g1,paste0(outfile,".graphml"),format = "graphml")
}