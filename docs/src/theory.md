# Fuzzy C-Regression

Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of heterogeneity. It extends the "Fuzzy C-Means" clustering algorithm to regression settings.

FCR can be used to estimate ["grouped fixed effects"](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA11319) (a constant term for each group) as well as heterogeneous coefficients.

In panel settings, each unit belongs to the same group over time, meaning that separate coefficients are estimated for each group-time pair.

## FCR objective function
Consider a linear model with grouped heterogeneity:

$$y = \sum_{g=1}^G\mu_{g(i)} \theta_{g}X+\varepsilon$$
where $X_{i}$ are covariates (which could simply be a constant), $\theta_{g}$ represent group-specific coefficients for groups $g=1,\ldots,G$, and $\mu_{g(i)}$ represent group weights for unit $i$. FCR is concerned with jointly estimating $\theta$ and $\mu$ for each group.

The FCR objective function takes the form:

$$L^{FCR}_m\left(\theta,\mu\right)=E\left[\sum_{g=1}^{G}\mu_{g}^{m}\left\Vert y-\theta_g x\right\Vert ^{2}\right]$$

where $m > 1$ is the regularization parameter. The weights are defined as

$$\mu_{g}\left(y,x;\theta,m\right)=\left(\sum_{h=1}^{G}\frac{\left\Vert y-\theta_g x\right\Vert ^{2/\left(m-1\right)}}{\left\Vert y-\theta_h x\right\Vert ^{2/\left(m-1\right)}}\right)^{-1} \text{ for } g=1,\ldots,G$$


Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors $||y-\theta_gX||$:

$$J^{FCR}\left(\theta\right)=E\left[\left(\sum_{g=1}^{G}\left\Vert y-\theta_g x\right\Vert ^{-2/\left(m-1\right)}\right)^{1-m}\right]$$

Thus, for fixed $m$, the FCR function is differentiable and can be written as a standard GMM problem. 

## Useful properties

  - __Fast:__ FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.

  - __Customizable:__ The “fuzziness” of the FCR groups is governed by the regularization parameter $m$, where group assignment becomes binary as $m \rightarrow 1$. Choosing $m$ allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.

  - __Inference:__ Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors. 
