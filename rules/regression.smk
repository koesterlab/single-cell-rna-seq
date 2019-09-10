
def get_regression_config(wildcards):
    return config["pairwise-regression"][wildcards.regression]

def get_constrain_celltypes(wildcards):
    return get_regression_config(wildcards).get("constrain-celltypes")


def get_gene_vs_gene_fits(wildcards):
    constrain_celltypes = get_constrain_celltypes(wildcards)
    constrained_markers = markers
    if constrain_celltypes:
        constrained_markers = markers.loc[markers["name"].isin(constrain_celltypes)]
    return expand("analysis/cellassign.{parent}.rds", parent=constrained_markers["parent"].unique())
    


rule gene_vs_gene:
    input:
        sce="analysis/normalized.batch-removed.rds",
        fits=get_gene_vs_gene_fits
    output:
        report("plots/gene-vs-gene/{gene_a}-vs-{gene_b}.{regression}.regression.pdf",
                   caption="../report/pairwise-regressions.rst",
                   category="Pairwise Regressions")
    params:
        min_gamma=config["celltype"]["min_gamma"],
        constrain_celltypes=get_constrain_celltypes,
        formula=lambda w: get_regression_config(w).get("formula", "y ~ x")
    log:
        "logs/gene-vs-gene/{gene_a}-vs-{gene_b}.{regression}.log"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-gene-gene-expression.R"
