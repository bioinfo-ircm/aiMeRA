% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netgraph.R
\name{netgraph}
\alias{netgraph}
\title{Network connectivity plot}
\usage{
netgraph(map,layout=igraph::layout_with_kk,pertu=NULL,inter=NULL,
      cutoff=NULL,main=NULL,digits=2,outfile=NULL,module.in=NULL,
      module.out=NULL, no.module=NULL,neg.col="red",pos.col="green", ...)
}
\arguments{
\item{map}{A list containing the connectivity map (link matrix) of a network, and local matrix responses
to perturbations, i.e., the output of ```mra```.}

\item{layout}{Either a data.frame containing the x and y coordinates and the color of the vertices in the network,
or a layout function according to the igraph package. Default uses the "layout_with_kk" igraph layout.}

\item{pertu}{The name of perturbations to be plotted as vertices in the network.}

\item{inter}{Confidence intervals calculated by ```interval```; connectivity coefficient with a confidence
interval that does not include 0 are deemed significant and marked by an asterisk.}

\item{cutoff}{Minimum value for a connectivity link coefficient between two modules for being plotted.}

\item{main}{A character string giving the title of the plot.}

\item{digits}{Number of digits to be plotted for each connectivity coefficient.}

\item{outfile}{Optional. A character string giving the name of the output file in graphML format.}

\item{module.in}{The name of one or more nodes for which only the incoming connectivity
from other nodes of the network should be plotted.}

\item{module.out}{The name of one or more nodes for which only the outcoming connectivity
to other nodes of the network should be plotted.}

\item{no.module}{The name of one or more nodes for which incoming and outcoming connectivity
should not be plotted.}

\item{neg.col}{A color to be used for arrows with negative connectivity coefficients}

\item{pos.col}{A color to be used for arrows with positive connectivity coefficients}

\item{...}{Arguments to be passed to igraph plot such as vertex and edges plotting parameters.}
}
\value{
If layout is an igraph function then a data.frame containing the coordinates of the layout selected is returned.
If outfile is not NULL, then the network is written in the graphML file format.
}
\description{
Plots a graphic representation of a biological network connectivity map calculated by MRA.
}
\examples{
data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
tec.av=list(data2sdmean(data[1:2])$mean,data2sdmean(data[3:4])$mean,data2sdmean(data[5:6])$mean)
data.mean=data2sdmean(tec.av)$mean
lb="E2"
data.rp=2*(data.mean[,colnames(data.mean)!=lb]-data.mean[,lb])/(data.mean[,colnames(data.mean)!=lb]+data.mean[,lb])
rules=c("E2+siLCoR->LCoR","E2+siRIP140->RIP140","Et->Luciferase","E2->0")
matp=read.rules(rules)
#The variance of each variable was estimated employing an estimator optimized for a
#small sample size from Statistical Process Control theory
#(Wheeler and Chambers, 1992; Harter, 1960). The standard deviation was computed for
#the global response matrices and stored into the sd.ex table which is included
#in the package)
inter=interval(data.rp,sd.tab=sd.ex,matp=matp,Rp=TRUE)
map=mra(data.rp,matp,Rp=TRUE,check=FALSE)
netgraph(map,inter=inter)
}
