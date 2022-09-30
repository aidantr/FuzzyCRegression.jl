# Fuzzy Clustering Regression

## FCR Objective 
Consider a linear model with grouped heterogeneity:

\begin{aligned}
y_{it} = \sum_{g=1}^G\mu_{g(i)} \theta_{g}X_i+\varepsilon_{it}
\end{aligned}
where $X_{it}$ are covariates (which could simply be a constant), $\theta_{g}$ represent group-specific coefficients for groups $g=1,\ldots,G$, and $\mu_{g(i)}$ represent group weights for unit $i$. Fuzzy Clustering Regression (FCR) is concerned with jointly estimating $\theta$ and $\mu$ for each group.

The FCR objective function takes the form:

\begin{aligned}
L^{FCR}_m\left(\theta,\mu\right)=E\left[\sum_{g=1}^{G}\mu_{g}^{m}\left\Vert y-\theta_g x\right\Vert ^{2}\right]
\end{aligned}

where $m > 1$ is the regularization parameter. The weights are defined as

\begin{aligned}
\mu_{g}\left(y,x;\theta,m\right)=\left(\sum_{h=1}^{G}\frac{\left\Vert y-\theta_g x\right\Vert ^{2/\left(m-1\right)}}{\left\Vert y-\theta_h x\right\Vert ^{2/\left(m-1\right)}}\right)^{-1},g=1,\ldots,G
\end{aligned}

Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors $||y=\theta_gX||$. Thus, for fixed $m$, the FCR function is differentiable and can be written as a standard GMM problem. This has two important implications:

## Useful Properties

1. FRC can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units into groups. While FCR remains quite efficient with large datasets, computation time for iterative approaches becomes prohibitive.

2. The “fuzziness” of the FCR groups is governed by the
regularization parameter $m$, where group assignment becomes binary as $m \rightarrow 1$. Choosing $m$ allows one to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it FCR can then recover the full distribution of effects.

3. Since FCR is a GMM problem, it's asymptotic properties follow from standard arguments and we are able to derive analytic standard errors. 
