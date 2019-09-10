Comparison of log counts for **{wildcards.gene_a}** with **{wildcards.gene_b}**. Each marker is a cell.
{% set constrain_celltypes = config["pairwise-regression"][wildcards.regression].get("constrain-celltypes", []) %}
{% if constrain_celltypes %}
Only cells determined as {{ constrain_celltypes|join(", ") }} are plotted.
{% endif %}
The blue line shows a performed linear regression with formula ``{{config["pairwise-regression"][wildcards.regression]["formula"]}}``. The determined formula, together with the coefficient of determination :math:`$R^2` (see `here <https://en.wikipedia.org/wiki/Coefficient_of_determination>`_) is shown within the plot.
