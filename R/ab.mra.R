#'@title Inference of the effect of a double perturbation on the modules
#'activity levels.
#'@description This function first computes the connectivity coefficients
#'between the modules of a biological network using standard MRA methodology.
#'Then, it infers the result of the simultaneous action of two of the
#'perturbations that were used for training the the MRA model.
#'@usage ab.mra(data, matp, pred=NULL, pert1, pert2, inval=c(-1,1), step=0.1,
#'Rp=FALSE, ab=TRUE)
#'@param data A data frame containing the experimental data in the required
#'            format for MRA computations (see \code{mra()}).
#'@param matp The perturbation matrix. Names of modules (rows) and
#'            perturbations (columns) must correspond to names of rows and
#'            columns in \code{tab}.
#'@param pred String. Name of the double perturbation to be inferred.
#'            This parameter should only be used when the perturbation to be
#'            inferred was actually measured experimentally, i.e., \code{pred}
#'            must appear as a column in \code{data}.
#'@param pert1 String. Name of the first elementary perturbation.
#'@param pert2 String. Name of the second elementary perturbation.
#'@param inval A two-value vector giving the lower and the upper limits of the
#'             intervals for the a and b coefficients (see Jimenez-Dominguez,
#'             et al., Sci Rep, 2021 for details).
#'@param step Number. Increment to cover \code{interval} for scanning the
#'            a and b coefficients values.
#' @param Rp Logical. TRUE if \code{tab} is the already computed global
#'           response matrix and not the original data. In the latter case,
#'           the predicted module activity will be expressed as its
#'           relative change compared to the basal condition. If
#'           \code{Rp=FALSE}, then the inferred values will be those of
#'           the module activity levels.
#'@param ab Logical. If TRUE, then the inferred values of a double perturbation
#'          are obtained by using the a and b coefficients. If FALSE, then
#'          a=1 and b=1 (see details).
#'@return List. If \code{pert} is provided and \code{ab=TRUE} then the inferred
#'value, the reference data, and the values for the a and b coefficients are
#'returned. If \code{pert} is provided and \code{ab=FALSE}, then the inferred
#'and the reference values are returned. If \code{pert} is not provided, then
#'only the inferred values are returned.
#'@details \code{ab.MRA()} inference relies on the assumption that the
#'MRA model is relevant and that the result of the combined perturbations
#'depends in an almost linear fashion from the corresponding individual
#'two perturbations.  See \code{Rb} parameter description
#'for the nature of the inferred activity levels (relative or absolute).
#'
#'Two coefficients (a and b) can be used to weigh the local responses to
#'each perturbation. If both are set to 1, it is assumed that the response
#'strengths remain identical to what they were under single perturbations.
#'Different values can be set if data are available to learn the relative
#'strengths under the combined action of the perturbations (see
#'Jimenez-Dominguez et al., Sci Rep, 2021 for more details and examples).
#'@export
#'@examples
#'# Inference of the effect of a double perturbation by two siRNAs (siRIP140
#'# and siLCoR) under the basal condition (E2-stimulated) in the
#'# ERa-RIP140-LCoR transcriptional network.
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'data.mean <- data2sdmean(data)$mean
#'rules <- c("Et->Luciferase", "E2+siRIP140->RIP140", "E2+siLCoR->LCoR",
#'           "E2->0")
#'matp <- read.rules(rules)
#'
#'# Inference without using the a and b coefficients
#'ab.mra(data.mean, matp=matp, pred="E2+siLCoR+siRIP140", pert1="E2+siLCoR",
#'       pert2="E2+siRIP140", ab=FALSE)
#'# Inference using the a and b coefficients
#'ab.mra(data.mean, matp=matp, pred="E2+siLCoR+siRIP140", pert1="E2+siLCoR",
#'        pert2="E2+siRIP140")
#'
#'# Inference without the experimental value of the reference and without
#'# the a and b coefficients
#'ab.mra(data.mean, matp=matp, pred="E2+siLCoR+siRIP140", pert1="E2+siLCoR",
#'       pert2="E2+siRIP140")
#'
#'# Inference of a biological module (GREB1), which does not have an individual
#'# perturbation and without using a,b
#'rules <- c("Et->Luciferase", "E2+siRIP140->RIP140", "E2+siLCoR->LCoR",
#'           "E2->0", "0->GREB1")
#'matp <- read.rules(rules)
#'ab.mra(data.mean, matp=matp, pred="E2+siLCoR+siRIP140", pert1="E2+siLCoR",
#'       pert2="E2+siRIP140", ab=FALSE)


ab.mra <- function(data, matp, pred=NULL, pert1, pert2, inval=c(-1,1), step=0.1,
                    Rp=FALSE, ab=TRUE) {
    if (!pert1 %in% colnames(matp) | !pert2 %in% colnames(matp))
        stop(paste0("pert1 and pert2 should match column ",
                    "names of the the perturbation matrix"))
    if (length(inval) != 2 | !is.numeric(inval))
        stop("The interval must be a vector with only two numerical values")
    res <- mra(data, matp, check=FALSE, Rp=Rp)
    rpa <- diag(res$local_matrix)
    coeff <- seq(inval[1], inval[2], step)
    lb <- colnames(matp)[colSums(matp) == 0]
    if (!is.null(pred)) {
        if (!pred %in% colnames(data))
            stop("pred should match a column name in data to be inferred")
        observed <- data[rownames(matp), pred]
    }

    aux <- colnames(res$local_matrix) %in% c(pert1, pert2)
    if (ab & !is.null(pred)) {
        glb <- data[rownames(matp), lb]
        c <- rep(0, ncol(res$local_matrix))
        errors <- vector()
        n <- length(coeff)
        for (k in coeff)
            for (l in coeff) {
                c[aux] <- c(k, l)
                infer <- solve(res$link_matrix, -rpa * c)
                if (!Rp)
                    infer <- (-glb * ((2 + infer) / (2 - infer)))
                errors <- c(errors, sqrt(sum((infer - observed)^2)))
            }
        err.mat <- matrix(errors, nrow=n, ncol=n, byrow=TRUE)
        ind.min <- arrayInd(which(err.mat == min(err.mat)), .dim=c(n, n))[1,]
        albeta <- coeff[ind.min]
        c[aux] <- albeta
        pred <- solve(res$link_matrix, -rpa * c)
        if (!Rp) {
            glb <- data[rownames(matp), lb]
            pred <- glb * (2 + pred) / (2 - pred)
        }
        all_coeff_total <- rep(0, ncol(res$local_matrix))
        all_coeff_total[aux] <- albeta
        names(all_coeff_total) <- rownames(res$link_matrix)
        list(coeff.ab=all_coeff_total, inferred=pred, observed=observed)
    }
    else {
        rpa[!aux] <- 0 # the weights in c (see above) are 0 or 1
        infer <- solve(res$link_matrix, -rpa)
        if (!Rp) {
            glb <- data[rownames(matp), lb]
            infer <- glb * (2 + infer) / (2 - infer)
        }
        if (!is.null(pred))
            list(inferred=infer, observed=observed)
        else
            list(inferred=infer)
    }
}
