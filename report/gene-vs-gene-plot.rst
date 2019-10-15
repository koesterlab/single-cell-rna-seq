Comparison of log counts for **{wildcards.gene_a}** with **{wildcards.gene_b}**. Each marker is a cell.
{% set regression = snakemake.params.regression %}
{% set correlation = snakemake.params.correlation %}
{% set constrain_celltypes = snakemake.params.constrain_celltypes %}
{% set dropout_threshold = snakemake.params.dropout_threshold %}
{% if constrain_celltypes %}
Only cells determined as {{ constrain_celltypes|join(", ") }} are plotted.
{% endif %}
{% if regression %}
The blue line shows a performed linear regression with formula ``{{config["pairwise-regression"][wildcards.regression]["formula"]}}``. The determined formula, together with the coefficient of determination :math:`$R^2` (see `here <https://en.wikipedia.org/wiki/Coefficient_of_determination>`_) is shown within the plot.
Only expressions above {{ dropout_threshold }} for both genes are considered.
{% endif %}
{% if correlation %}
The {{ correlation }} correlation coefficient together with the corresponding p-value is shown within the plot, along with a linear regression (blue line) including the 95% confidence interval (blueish area).
Only expressions above {{ dropout_threshold }} for both genes are considered.
{% endif %}
All expressions below {{ dropout_threshold }} for any of the two genes (shown in grey) are considered as dropouts.
