#' @title Create a perturbation matrix from a set of perturbation rules.
#' @description A function that takes a set of perturbation rules as input and
#'              creates a perturbation matrix to be used for MRA analysis.
#' @usage read.rules(obj)
#' @param obj Either the path to a perturbation rules file without header,
#'            or a vector of strings specifying all the perturbation rules
#'            (see details).
#' @details Perturbations rules are a set of strings defining the action of
#'          perturbations upon modules in a network. In both cases, modules
#'          and perturbations must have the same names as the names of rows and
#'          columns in the experimental data tables (case sensitive).
#'
#'The rule syntax is  "Perturbation->Module" for a single
#'perturbation, and "Perturbation->0" to specify the basal condition.
#'"0->Module" can be used to specify modules that were not perturbed,
#'but for which connectivity coming from other modules should be retrieved.
#'At least two rules and basal condition must be provided.
#'
#' @export
#' @importFrom data.table fread
#' @importFrom data.table is.data.table
#' @importFrom data.table setDF
#' @examples
#'rules <- c("Et->Luciferase", "E2+siRIP140->RIP140", "E2+siLCoR->LCoR",
#'           "E2->0")
#'read.rules(rules)

read.rules <- function(obj) {
    if (is.character(obj) & length(obj) == 1)
        rules <- data.table::fread(obj, stringsAsFactors=FALSE,
                                   data.table=FALSE, header=FALSE)$V1
    else if (is.data.frame(obj)) {
        if (ncol(obj) > 1)
            stop("Only a 1-column table was expected")
        if(is.data.table(obj))
            obj <- setDF(obj)
        rules <- obj[,1]
    }
    else if (is.vector(obj))
        rules <- obj
    else
        stop("Incorrect input format")

    if (length(rules) < 2)
        stop("Not enough perturbations for the MRA analysis")
    if (!all(grepl("->", rules)))
        stop(paste("Incorrect syntax in : ",paste(rules[!grepl("->",rules)],
                                                  collapse=", "),sep=""))
    if (!any(grepl("*->0$",rules)))
        stop("The basal condition is not specified")
    if (sum(grepl("*->0$",rules)) > 1)
        stop("The basal condition is defined multiple times")

    if (any(grepl(" ", rules)))
        rules <- gsub(" ", "", rules)
    if (any(grepl("^0->*", rules)))
        bid <- grepl("^0->*", rules)
    else
        bid <- NULL
    rules <- strsplit(rules, "->")
    rules.from <- vapply(rules, function(x) x[1], FUN.VALUE="")
    rules.to <- vapply(rules, function(x) x[2], FUN.VALUE="")
    if (any(duplicated(rules.from[rules.from != "0"])))
        stop("At least one perturbation affects more than one module")
    if (any(duplicated(rules.to)))
        stop("At least one module is affected by more than one perturbation")

    if (!is.null(bid)) {
        mat <- matrix(0, nrow=length(rules) - (1 + sum(bid)),
                 ncol=length(rules) - (1 + sum(bid)),
                 dimnames=list(rules.to[rules.to != "0" & rules.from != "0"],
                               rules.from[rules.to != "0" & rules.from != "0"]))
        aux <- matrix(0, nrow=sum(bid), ncol=ncol(mat))
        rownames(aux) <- rules.to[rules.from == "0"]
        mat <- rbind(mat, aux)
    }
    else
        mat <- matrix(0, nrow=length(rules) - 1, ncol=length(rules) - 1,
                      dimnames=list(rules.to[rules.to != "0"],
                                    rules.from[rules.to != "0"]))
    diag(mat) <- 1
    mat <- cbind(mat, 0)
    colnames(mat)[ncol(mat)] <- rules.from[rules.to == "0"]
    mat
}
