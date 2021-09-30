#'@title Standard deviation and average computation for a data set.
#'@description Compute the standard deviation and average values of a data set.
#'Namely, if n replicates were generated,
#'then n data tables should be provided to \code{data.setup()} and
#'\code{data2sdmean()} computes the average and standard deviation
#'of each module over the replicates of each condition (or perturbation).
#'@usage data2sdmean(x)
#'@param x A list containing 2 or more tables.
#'@return A list containing two matrices, the standard deviations and mean
#'        values for each variable in x.
#'@export
#'@examples
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'data2sdmean(data)


data2sdmean <- function(x) {
    if(!is.list(x))
        stop("Input is not a list of data tables")
    m <- Reduce("+", x) / length(x)
    sd <- sqrt(Reduce("+", lapply(x, function(y) (y - m)^2)) / (length(x) - 1))
    return(list(sd=as.matrix(sd), mean=as.matrix(m)))
}
