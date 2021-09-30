#' @title Confidence interval estimation for MRA connectivity coefficients
#' @description A function that estimates confidence intervals of MRA
#'              connectivity coefficients based on replicates. It uses
#'              bootstrapping.
#' @usage interval(tab, mean=0, sd.tab, matp, n=10000, digits=2, alpha=0.05)
#' @param tab The global response matrix.
#'            by \code{data.setup()}.
#' @param mean mean of the normal distribution for creating the noisy matrices.
#' @param sd.tab A data table containing the standard deviation values of
#'               replicates for sampling global response matrices in the
#'               bootstrap.
#' @param matp A perturbation rule table, which rows correspond to MRA modules
#'             and columns to perturbations.
#' @param n Number of samples for the bootstrap.
#' @param digits If different from -1 the interval boundaries are rounded at
#'               this number of digits.
#' @param alpha significance of the confidence interval.
#' @return A list containing the upper and lower values of the confidence
#'         interval of each connectivity coefficient.
#' @export
#' @importFrom stats rnorm
#' @examples
#'# Confidence intervals are estimated based on the variability in the global
#'# responses (R). Basal condition is defined as estradiol (E2) stimulation
#'# for the ER-RAR-LCoR-RIP140 example provided with the package.
#'# We first average over the two technical replicates of each biological
#'# replicate before averaging on the biological replicates themselves.
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'tec.av <- list(data2sdmean(data[1:2])$mean, data2sdmean(data[3:4])$mean,
#'               data2sdmean(data[5:6])$mean)
#'data.mean <- data2sdmean(tec.av)$mean
#'data.rp <- global.matrix(data.mean, "E2")
#'rules <- c("E2+siLCoR->LCoR", "E2+siRIP140->RIP140", "Et->Luciferase",
#'           "E2->0")
#'matp <- read.rules(rules)
#'# Estimation using the inter-transfection variability
#'intertransf.sd <- data2sdmean(global.matrix(tec.av, lb="E2"))$sd / sqrt(3)
#'interval(data.rp, sd.tab=intertransf.sd, matp=matp)
#'# The above result can be replaced by estimates obtained using the intra-
#'# transfection variability, which might be arguably more appropriate in this
#'# example data set (see Jimenez-Dominguez et al., Sci Rep,
#'# 2021 for more explanations).
#'# This illustrates how aiMeRA users can substitute their own estimates of
#'# data variability. In the particular case, the standard deviation of each
#'# variable was estimated employing an estimator optimized
#'# for small sample sizes from Statistical Process Control (SPC) theory
#'# (Wheeler and Chambers, 1992; Harter, 1960).
#'Rbs <- global.matrix(data, lb="E2")
#'d2nr3 <- 1.128
#'intratransf.sd <- (abs(Rbs[[1]] - Rbs[[2]]) +
#'                   abs(Rbs[[3]] - Rbs[[4]]) +
#'                   abs(Rbs[[5]] - Rbs[[6]])
#'                  ) / (3 * d2nr3)
#'interval(data.rp, sd.tab=intratransf.sd, matp=matp)

interval <- function(tab, mean=0, sd.tab, matp, n=10000, digits=2, alpha=0.05) {
    if(!is.matrix(sd.tab))
        sd.tab <- as.matrix(sd.tab)
    lb <- colnames(matp)[colSums(matp) == 0]
    pt <- colnames(matp)[colnames(matp) != lb]
    genes <- rownames(matp)
    tab <- tab[genes, pt]
    sd.tab <- sd.tab[genes,pt]
    rlist <- list()
    rplist <- list()
    for (i in seq_len(n)) {
        xnoise <- tab + matrix(vapply(as.vector(sd.tab),
                               function(x)rnorm(1, mean, x),
                               FUN.VALUE=1.0
                               ), ncol=ncol(sd.tab))
        res <- mra(xnoise, matp, check=FALSE, Rp=TRUE)
        diag(res$link_matrix) <- NA
        r <- as.vector(res$link_matrix)
        rp <- diag(res$local_matrix)
        r <- r[!is.na(r)]
        if (i == 1) {
            rlist <- as.list(r)
            rplist <- as.list(rp)
        }
        else {
            rlist <- lapply(seq_len(length(r)), function(j) c(rlist[[j]], r[j]))
            rplist <- lapply(seq_len(length(rp)),
                            function(j) c(rplist[[j]], rp[j]))
        }
    }
    noms <- expand.grid(rep(list(genes), 2))
    noms <- noms[noms[[1]] != noms[[2]], 2:1]
    names(rlist) <- apply(noms, 1, function(x) paste0(x, collapse = "->"))
    names(rplist) <- paste0(colnames(res$local_matrix), "->",
                            rownames(res$local_matrix))
    inter <- c(lapply(rlist, function(x) c(stats::quantile(x,prob=0.5*alpha),
                                       stats::quantile(x,prob=1-0.5*alpha))),
               lapply(rplist, function(x) c(stats::quantile(x,prob=0.5*alpha),
                                       stats::quantile(x,prob=1-0.5*alpha))))
    if (digits > -1)
        inter <- lapply(inter, function(x) round(x, digits=digits))
    inter
}
