#'@title Network connectivity plot
#
#'@description Plot a graphical representation of a biological network
#'             connectivity as computed by MRA.
#'
#'@usage netgraph(map, layout=igraph::layout_with_kk, pertu=NULL, inter=NULL,
#'       cutoff=NULL, main=NULL, digits=2, outfile=NULL, module.in=NULL,
#'       module.out=NULL, no.module=NULL, neg.col="red", pos.col="green", ...)
#'@param map A list containing the topology definition of a network, and
#'           MRA local response matrix, i.e., the output of \code{mra}.
#'@param layout Either a data.frame containing the x and y coordinates and the
#'              colors of the network vertices, or a layout function according
#'              the igraph package. Default layout is provided by the
#'              \code{layout_with_kk()} function of igraph.
#'@param pertu The name of the perturbations to be plotted as vertices in the
#'             network.
#'@param inter Confidence intervals estimated by \code{interval()}; connectivity
#'             coefficients featuring confidence intervals that do not include 0
#'             are deemed significant and marked by an asterisk.
#'@param cutoff Minimum value for a connectivity coefficient to be plotted.
#'@param main A character string giving the title of the plot.
#'@param digits Number of digits to be plotted for each connectivity
#'              coefficient.
#'@param module.in The name of one or more modules for which only the incoming
#'              connectivity from the other modules should be plotted.
#'@param module.out The name of one or more modules for which only the outcoming
#'              connectivity to the other modules should be plotted.
#'@param no.module The name of one or more modules to which the plot should be
#'              restricted.
#'@param neg.col Negative connectivity arrow color.
#'@param pos.col Positive connectivity arrow color.
#'@param outfile Optional. An output file name to write in graphML format.
#'@param ... Arguments to be passed to igraph plot such as vertex and edge
#'           plotting parameters.
#'@return If layout is an igraph function, then a data.frame containing the
#'        coordinates of the selected layout is returned.
#'        Provided \code{outfile} is not \code{NULL}, the network is written in
#'        graphML file format.
#'@export
#'@importFrom igraph V
#'@importFrom igraph graph
#'@importFrom igraph set_edge_attr
#'@importFrom igraph write.graph
#'@importFrom igraph %>%
#'@importFrom graphics plot
#'@examples
#'data <- data.setup(list(estr1_A, estr1_B, estr2_A, estr2_B, estr3_A, estr3_B))
#'tec.av <- list(data2sdmean(data[1:2])$mean, data2sdmean(data[3:4])$mean,
#'               data2sdmean(data[5:6])$mean)
#'data.mean <- data2sdmean(tec.av)$mean
#'lb <- "E2"
#'data.rp <- 2 * (data.mean[,colnames(data.mean) != lb] - data.mean[,lb]) /
#'                   (data.mean[,colnames(data.mean) != lb] + data.mean[,lb])
#'rules <- c("E2+siLCoR->LCoR", "E2+siRIP140->RIP140", "Et->Luciferase",
#'           "E2->0")
#'matp <- read.rules(rules)
#'# The variance of each variable was estimated employing an estimator optimized
#'# for small sample sizes from Statistical Process Control theory
#'# (Wheeler and Chambers, 1992; Harter, 1960). The standard deviation was
#'# computed for the global response matrix coefficients and stored into the
#'# sd.ex table which is part of this package.
#'Rbs <- global.matrix(data, lb="E2")
#'d2nr3 <- 1.128
#'intratransf.sd <- (abs(Rbs[[1]] - Rbs[[2]]) +
#'                   abs(Rbs[[3]] - Rbs[[4]]) +
#'                   abs(Rbs[[5]] - Rbs[[6]])
#'                  ) / (3 * d2nr3)
#'inter <- interval(data.rp, sd.tab=intratransf.sd, matp=matp)
#'map <- mra(data.rp, matp, Rp=TRUE, check=FALSE)
#'netgraph(map, inter=inter)

netgraph <- function(map, layout=igraph::layout_with_kk, pertu=NULL, inter=NULL,
                    cutoff=NULL, main=NULL, digits=2, outfile=NULL,
                    module.in=NULL, module.out=NULL, no.module=NULL,
                    neg.col="red", pos.col="green", ...) {
    map <- lapply(map, round, digits=digits)
    coeff <- numeric()
    noms <- expand.grid(list(rownames(map$link_matrix),
                            colnames(map$link_matrix)))
    noms <- noms[!apply(noms, 1, function(x)any(duplicated(x))),
                c("Var2","Var1")]
    if (!is.null(module.in) & is.null(module.out))
        noms <- noms[noms$Var1 %in% module.in,]
    else if (is.null(module.in) & !is.null(module.out))
        noms <- noms[noms$Var2 %in% module.out,]
    else if (!is.null(module.in) & !is.null(module.out))
        noms <- noms[noms$Var2 %in% module.out | noms$Var1 %in% module.in,]
    if (!is.null(no.module))
        noms <- noms[!(noms$Var2 %in% no.module|noms$Var1 %in% no.module),]
    for (i in seq_len(nrow(noms)))
        coeff <- c(coeff, map$link_matrix[noms[i, 2], noms[i, 1]])
    if (!is.null(cutoff)) {
        if(cutoff < max(abs(coeff)))
            coeff[abs(coeff) < cutoff] <- 0.0
        else
            stop("Cutoff too large. At least one network connnexion must exist")
    }
    if (any(coeff == 0.0)) {
        bid <- apply(map$global_matrix, 2, function(x) which(x == 0.0))
        bid <- unlist(lapply(bid,function(x) length(x)))
        noms <- noms[coeff != 0.0,]
        coeff <- coeff[coeff != 0.0]
    }
    names(coeff) <- apply(noms, 1, function(x) paste0(x, collapse = "->"))
    noms <- as.vector(t(noms))
    if (!is.null(pertu)) {
        if (!all(pertu %in% colnames(map$local_matrix)))
            stop("Perturbations to plot are not in the map object")
        mols <- rownames(map$local_matrix)[
                            colnames(map$local_matrix) %in% pertu]
        loc <- diag(map$local_matrix)[colnames(map$local_matrix) %in% pertu]
        names(loc) <- paste0(pertu, "->", mols)
        noms <- c(noms, rbind(pertu, mols))
        coeff <- c(coeff, loc)
    }
    width <- log(abs(coeff / min(abs(coeff[coeff != 0]))))
    width <- log(abs(coeff / 0.01))
    ncols <- ifelse(coeff < 0, neg.col, pos.col)
    g1 <- igraph::graph(noms, directed=TRUE)
    g1 <- g1 %>% set_edge_attr(name="coeff", value=coeff)
    if (!is.null(inter)) {
        bin <- lapply(inter[!unlist(lapply(inter, function(x)
                        sum(round(x, digits=digits)) == 0))],
                        round, digits=digits)
        bin <- !vapply(bin, function(x)
                                (x[1] < 0 & x[2] > 0) | (x[1] > 0 & x[2] < 0),
                        FUN.VALUE=TRUE)
        names(bin) <- gsub(".5%", "", names(bin))
        g1 <- g1 %>% set_edge_attr(name="sig",
                            value=ifelse(bin, "*", "")[seq_len(length(coeff))])
        coeff <- round(coeff, digits=digits)
        coeff[names(bin)][bin] <- paste0(coeff[names(bin)][bin], "*")
    }
    if (is.function(layout)) {
        coords <- do.call(layout,list(g1))
        plot(g1, edge.label=coeff, layout=coords, ... , vertex.size=10,
             edge.width=width, edge.curved=0.2, edge.color=ncols,
             edge.arrow.size=0.5, main=main)
        rownames(coords) <- igraph::V(g1)$name
        colnames(coords) <- c("X", "Y")
        if(!is.null(outfile))
            igraph::write.graph(g1, paste0(outfile, ".graphml"),
                                format = "graphml")
        return(invisible(coords))
    }
    plot(g1, edge.label=coeff, ... , vertex.size=10, edge.width=width,
         edge.curved=0.2, edge.color=ncols,
         edge.arrow.size=1, layout=layout[igraph::V(g1)$name,
                                                    c("X", "Y")]*0.1,
         rescale=FALSE, edge.label.cex=1, vertex.label.cex=1,
         vertex.color=layout[igraph::V(g1)$name,"Col"], main=main)
    if (!is.null(outfile))
        igraph::write.graph(g1, paste0(outfile, ".graphml"), format = "graphml")
}
