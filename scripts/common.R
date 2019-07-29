assign_celltypes <- function(cellassign_fit, sce, min_gamma) {
    max_gamma <- apply(cellassign_fit$mle_params$gamma, 1, max)
    # only consider cases where we are certain
    cell_type <- cellassign_fit$cell_type[max_gamma >= min_gamma,, drop=FALSE]
    # assign determined cell types
    colData(sce)[rownames(cell_type), "celltype"] <- sapply(cell_type, as.character)

    sce
}
