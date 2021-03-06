% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interval.R
\name{interval}
\alias{interval}
\title{Confidence interval calculation for MRA connectivity coefficients}
\usage{
interval(tab,mean=0,sd.tab,matp,Rp=FALSE,n=10000)
}
\arguments{
\item{tab}{A data table containing the experimental data in the format output by ```data.setup```.}

\item{mean}{mean of the normal distribution for creating the noisy matrices.}

\item{sd.tab}{A data table containing the standard deviation values of replicates for creating the noisy matrices.}

\item{matp}{A perturbation rule table, with rows corresponding to the MRA modules and columns to perturbations.}

\item{Rp}{Logical. TRUE if ```tab``` is the calcuated global response matrix.}

\item{n}{Number of samples for the bootstrap algorithm.}
}
\value{
A list containing the upper and lower values of the confidence interval of each connectivity coefficient.
}
\description{
A function that estimates confidence intervals of MRA connectivity coefficients
based on replicates. It uses a bootstrap algorithm.
}
\examples{
#Confidence intervals are obtained from the variability within the global
#reponse matrices (R). The basal condition is estradiol (E2)
#'#We first average the technical replicates of each biological replicate
#and then only keep the (averaged) biological replicates for the calculation
data=data.setup(list(estr1_A,estr1_B,estr2_A,estr2_B,estr3_A,estr3_B))
tec.av=list(data2sdmean(data[1:2])$mean,data2sdmean(data[3:4])$mean,data2sdmean(data[5:6])$mean)
data.mean=data2sdmean(tec.av)$mean
data.rp=global.matrix(data.mean,"E2")
rules=c("E2+siLCoR->LCoR","E2+siRIP140->RIP140","Et->Luciferase","E2->0")
matp=read.rules(rules)
#The variance of each variable was estimated employing an estimator optimized for a
#small sample size from Statistical Process Control theory
#(Wheeler and Chambers, 1992; Harter, 1960). The standard deviation was computed for
#the global response matrices and stored into the sd.ex table which is included
#in the package)
interval(data.rp,sd.tab=sd.ex,matp=matp,Rp=TRUE)
}
