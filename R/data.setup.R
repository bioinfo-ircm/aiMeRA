#'@title Data preparation for Modular Response Analysis
#'@description Check for possible errors in tables of experimental data and in
#'             the perturbation matrix. Return data properly formatted for the
#'             \code{mra()} function.
#'@param obj A character string specifying the path to a folder that contains
#'           the experimental data tables and the perturbation rules,
#'           or a list containing several data tables, or a single data table.
#'@param pert.tab Optional. The perturbation rule table if not provided in
#'                \code{obj}.
#'@param dec The decimal separator as in \code{base::read.csv} (default ".").
#'@return If more than one data table were checked, then a list of matrices in
#'        the correct format for \code{mra()} is returned. If a single table
#'        was checked, then a single matrix in the correct format is returned.
#'        \code{data.setup()} also checks the perturbation rule table.
#'@export
#'@importFrom data.table fread
#'@examples
#'# Example using the q-PCR data of Jimenez-Dominguez et al., Sci Rep, 2021
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'#Examples using a gene network from Kholodenko et.al., PNAS, 2002
#'data.A <- data.setup(gene.network.A)
#'data.B <- data.setup(gene.network.B)

data.setup <- function(obj, pert.tab=NULL, dec=".") {
    if (!is.list(obj) & is.vector(obj)) {
        if (length(obj) > 1)
            stop("Only one path to the data and rule tables is allowed")
        temp <- list.files(obj)
        if (length(temp) < 2)
            stop(paste0("At least one data table and a perturbation rule file ",
                        "are required"))
        data <- lapply(paste0(obj,temp), data.table::fread, data.table=FALSE,
                        dec=dec)
        if (any(vapply(data, ncol, FUN.VALUE=1) == 1)) {
            if (sum(vapply(data, ncol, FUN.VALUE=1) == 1) > 1)
                stop("Only one rule table is allowed")
            if (all(vapply(unlist(data[vapply(data, ncol, FUN.VALUE=1) == 1]),
                            is.character, FUN.VALUE=TRUE))) {
                pert.tab <- c(colnames(data[vapply(data, ncol,
                                                FUN.VALUE=1) == 1][[1]]),
                            unlist(data[vapply(data, ncol, FUN.VALUE=1) == 1]))
                data <- data[vapply(data, ncol, FUN.VALUE=1) != 1]
            }
            else
                stop("Wrong format in the file with perturbation rules")
        }
        else
            stop("File with perturbation rules was not found")
    }
    else if (is.data.frame(obj) | is.matrix(obj) | is.data.table(obj))
        data <- list(obj)
    else
        data <- obj
    if (any(vapply(data, ncol, FUN.VALUE=1) < 4))
        stop(paste0("Columns in data tables must contain at least two ",
                    "perturbations for two modules, modules names, ",
                    "and the basal condition"))
    if (!any(vapply(data, function(x) is.character(x[, 1]) |
                    is.character(colnames(x)), FUN.VALUE=TRUE)))
        stop("First column of each data table must contain the modules name")

    if (all(vapply(data, function(x) is.character(x[, 1]), FUN.VALUE=TRUE)))
        data <- lapply(data, function(x) {row.names(x) <- as.character(x[, 1]);
                                            x[, -1]})
    if (!is.null(pert.tab)) {
        if ((is.data.frame(pert.tab) | is.matrix(pert.tab))) {
            temp <- vapply(colnames(pert.tab), function(x)
                paste(x, rownames(pert.tab)[!!pert.tab[, x]], sep="->"),
                FUN.VALUE="")
            temp[colSums(pert.tab) == 0] <- paste0(temp[colSums(pert.tab) == 0],
                                                    "0")
            pert.tab <- temp
        }
        if (any(grepl(" ", pert.tab)))
            pert.tab <- gsub(" ", "", pert.tab)
        aux <- strsplit(pert.tab, "->")
        aux.from <- vapply(aux, function(x) x[1], FUN.VALUE="")
        aux.to <- vapply(aux, function(x) x[2], FUN.VALUE="")
        if (!any(grepl("*->0", pert.tab)))
            stop("The basal condition is not specified")
        if (sum(grepl("*->0$", pert.tab)) > 1)
            stop("The basal condition is defined more than once")
        if (any(duplicated(aux.from[aux.from != "0"])))
            stop("At least one perturbation affects more than one module")
        if (any(duplicated(aux.to)))
            stop(paste0("One or more modules are affected by more than ",
                        "one perturbation"))
        indices <- !grepl("^0->*", pert.tab)
        if (!all(vapply(aux[indices],
                        function(x) x[1] %in% colnames(data[[1]]),
                        FUN.VALUE=logical(sum(indices)))))
            stop("One perturbation name does not match names in a data table")
        indices <- !grepl("*->0$", pert.tab)
        if (!all(vapply(aux[indices],
                  function(x) x[2] %in% switch(2 - is.character(data[[1]][, 1]),
                                        data[[1]][, 1], rownames(data[[1]])),
                  FUN.VALUE=logical(sum(indices))
                 )))
            stop(paste0("One module name in the perturbation rule file does ",
                        "not match modules names in the data tables"))
    }
    if (length(data) > 1) {
        aux <- Reduce(paste, lapply(data, colnames))
        aux <- strsplit(aux, " ")
        if(!all(vapply(aux, function(x) all(x[1] == x), FUN.VALUE=TRUE)))
            stop(paste0("At least the name of one column is not the same ",
                        "between data tables"))
        aux <- Reduce(paste, lapply(data, rownames))
        aux <- strsplit(aux, " ")
        if(!all(vapply(aux, function(x) all(x[1] == x), FUN.VALUE=TRUE)))
            stop(paste0("At least the name of one row is not the same between ",
                        "data tables"))
        data <- lapply(data, as.matrix)
        return(data)
    }
    as.matrix(data[[1]])
}
