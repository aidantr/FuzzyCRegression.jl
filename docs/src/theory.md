# Fuzzy Clustering Regression

Consider the linear model:
$$ y_{it} = x_{it}'\theta + \sum^G\alpha_{g_it}+\varepsilon_{it}$$
where $x_{it}$ are covariates uncorrelated with the error term and $\alpha_{g_it}$ represent group-specific unobservables for groups $g=1,\ldots,G$. Fuzzy Clustering Regression (henceforth FCR) is concerned with estimating $(\theta,\alpha_g)$.

The FCR objective function takes the form:
