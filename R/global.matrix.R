#'@title Global response matrix computation
#'
#'@description Compute the matrix containing the relative observed changes upon
#'             perturbation.
#'
#'@usage global.matrix(data, lb)
#'@param data A data.frame containing experimental data in a specific format.
#'@param lb Basal condition to compute the relative changes in module activity.
#'@return A data.frame containing the relative activity changes.
#'@export
#'@examples
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'Rp <- global.matrix(data, lb="E2")


global.matrix <- function(data, lb) {
    if (is.list(data))
        lapply(data, function(x)
            2 * (x[,colnames(x) != lb] - x[, lb]) /
                            (x[, colnames(x) != lb] + x[, lb]))
    else
        2 * (data[, colnames(data) != lb] - data[, lb]) /
                            (data[, colnames(data) != lb] + data[, lb])
}
