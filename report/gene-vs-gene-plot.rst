Comparison of log counts for **{wildcards.gene_a}** with **{wildcards.gene_b}**. Each marker is a cell.
{% set regression = snakemake.params.regression %}
{% set correlation = snakemake.params.correlation %}
{% set constrain_celltypes = snakemake.params.constrain_celltypes %}
{% if constrain_celltypes %}
Only cells determined as {{ constrain_celltypes|join(", ") }} are plotted.
{% endif %}
{% if regression %}
The blue line shows a performed linear regression with formula ``{{config["pairwise-regression"][wildcards.regression]["formula"]}}``. The determined formula, together with the coefficient of determination :math:`$R^2` (see `here <https://en.wikipedia.org/wiki/Coefficient_of_determination>`_) is shown within the plot.
{% endif %}
{% if correlation %}
The {{ correlation }} correlation coefficient together with the corresponding p-value is shown within the plot, along with a linear regression (blue line) including the 95% confidence interval (blueish area).
{% endif %}
