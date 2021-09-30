#' @title Modular Response Analysis model inference function
#'
#' @description Calculation of MRA connectivity coefficients.
#' @usage mra(tab, matp, check=TRUE, Rp=FALSE)
#' @param tab A data.frame containing experimental data in a specific format
#'            (see details).
#' @param matp The perturbation matrix. Names of modules (rows) and
#'             perturbations (columns) must correspond to names of rows and
#'             columns in tab.
#' @param check Logical. Should the dataset and perturbation matrix be checked
#'              for input errors?
#' @param Rp Logical. TRUE if \code{tab} is the already computed global
#'           response matrix and not the original data.
#' @return A list containing the connectivity map, the local responses matrix,
#'         the network responses matrix to perturbations,
#'         and the basal values for all the modules.
#' @details It is assumed that one perturbation only affects one module of the
#'          network. This is coded as binary values in the perturbation matrix
#'          \code{matp}. It is further assumed that row names in \code{tab}
#'          are the module names, and column names are the perturbation names.
#' @export
#' @examples
#'# Create the connectivity map between 2 transcriptional nuclear coregulators
#'# (RIP140 and LCoR) and estrogen receptor alpha transcriptional activity as
#'# reported by a luciferase gene. q-PCR data are stored in the package files
#'# used below calling the function data.setup.
#'# The model is obtained under E2 stimulation as basal condition.
#'
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'sd.mean <- data2sdmean(data)
#'rules <- c("Et->Luciferase", "E2+siRIP140->RIP140",
#'           "E2+siLCoR->LCoR", "E2->0")
#'matp <- read.rules(rules)
#'mra(sd.mean$mean, matp)

mra <- function (tab, matp, check=TRUE, Rp=FALSE) {
    if (check)
        data.setup(tab, matp)
    if (!is.matrix(matp))
        matp <- as.matrix(matp)
    lb <- colnames(matp)[colSums(matp) == 0]
    if(!Rp) {
        glb <- tab[rownames(matp), lb]
        gl <- tab[rownames(matp), colnames(matp)[colnames(matp) != lb]]
        Rp <- as.matrix(2 * (gl - glb) / (gl + glb))
    }
    else
        Rp <- as.matrix(tab[rownames(matp),
                            colnames(matp)[colnames(matp) != lb]])
    rownames(Rp) <- rownames(matp)
    colnames(Rp) <- colnames(matp)[colSums(matp) != 0]
    if (any(rowSums(matp) == 0)) {
        aux <- matrix(0, ncol=sum(rowSums(matp) == 0), nrow=nrow(Rp))
        colnames(aux) <- rep("bid", sum(rowSums(matp) == 0))
        rownames(aux) <- rownames(Rp)
        aux[rowSums(matp) == 0,] <- diag(ncol(aux))
        Rp <- cbind(Rp, aux)
        matp <- cbind(matp, aux)
    }
    matp <- matp[,colSums(matp) != 0]
    if (!all(matp == diag(rep(1,ncol(matp)))))
        stop(paste0("Problem in perturbation rules: matp first columns should ",
             "be the identity matrix"))
    Rpinv <- solve(Rp)
    rp <- diag(1 / diag(Rpinv))
    rownames(rp) <- rownames(Rp)
    colnames(rp) <- colnames(matp)
    r <- -rp %*% Rpinv
    list(link_matrix=r, global_matrix=Rp, local_matrix=rp)
}
