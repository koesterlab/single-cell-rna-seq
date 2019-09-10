# Copyright 2018 Johannes Köster.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(scater)
library(scran)
library(ggsci)
library(ggpubr)
source(file.path(snakemake@scriptdir, "common.R"))

gene_a <- snakemake@wildcards[["gene_a"]]
gene_b <- snakemake@wildcards[["gene_b"]]

aes_name <- function(name) {
    paste0("`", name, "`")
}

sce <- readRDS(snakemake@input[["sce"]])

for(cellassign_fit in snakemake@input[["fits"]]) {
    cellassign_fit <- readRDS(cellassign_fit)
    sce <- assign_celltypes(cellassign_fit, sce, snakemake@params[["min_gamma"]])
}

# handle constrain celltypes
constrain_celltypes <- snakemake@params[["constrain_celltypes"]]
if(!is.null(constrain_celltypes)) {
    celltypes <- constrain_celltypes
    sce <- sce[, colData(sce)$celltype %in% celltypes]
}

# plot t-SNE
pdf(file=snakemake@output[[1]], width = 4, height = 4)
data <- as.data.frame(t(logcounts(sce[c(gene_a, gene_b), ])))
data <- data[apply(data, 1, min) >= 0, ]
formula <- as.formula(snakemake@params[["formula"]])
ggplot(data, aes_string(x=aes_name(gene_a), y=aes_name(gene_b))) +
    geom_point(shape=1) +
    geom_smooth(method="lm", formula = formula) +
    stat_regline_equation(label.x.npc = "left", label.y.npc = "top", formula = formula,
			  aes(label = paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
    theme_classic()

dev.off()
