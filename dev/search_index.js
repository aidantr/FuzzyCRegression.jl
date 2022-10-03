var documenterSearchIndex = {"docs":
[{"location":"theory/#Fuzzy-C-Regression","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"","category":"section"},{"location":"theory/#FCR-Objective","page":"Fuzzy C-Regression","title":"FCR Objective","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"y_it = sum_g=1^Gmu_g(i) theta_gtX_it+varepsilon_it","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where X_it are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g(i) represent group weights for unit i. Fuzzy Clustering Regression (FCR) is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"The FCR objective function takes the form:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g xrightVert ^2right","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"mu_gleft(yxthetamright)=left(sum_h=1^GfracleftVert y-theta_g xrightVert ^2left(m-1right)leftVert y-theta_h xrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"J^FCRleft(thetaright)=Eleftleft(sum_g=1^GleftVert y-theta_g xrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"theory/#Useful-Properties","page":"Fuzzy C-Regression","title":"Useful Properties","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.\nThe “fuzziness” of the FCR groups is governed by the regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows one to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.\nSince FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors. ","category":"page"},{"location":"package/#FuzzyCRegression.jl-Tutorial","page":"Julia Implementation","title":"FuzzyCRegression.jl Tutorial","text":"","category":"section"},{"location":"package/#Installation","page":"Julia Implementation","title":"Installation","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"Pkg.add(\"FuzzyCRegression\")","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"will install the package and its dependencies.","category":"page"},{"location":"package/#Fitting-the-FCR-model","page":"Julia Implementation","title":"Fitting the FCR model","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"An FCR model is fit using fit(y,X,G,m,...). The arguments are","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"y a vector holding values of the dependent variable\nX a matrix holding values of the independent variable(s) that will have heterogeneous coefficients (the default is a constant term, which will estimate a fixed effect for each group)\nZ a matrix holding values of the independent variable(s) that will have homogeneous coefficient\nG number of groups\nm regularization parameter (greater than 1), where group assignment becomes binary as m rightarrow 1","category":"page"},{"location":"package/#Methods-applied-to-fitted-models","page":"Julia Implementation","title":"Methods applied to fitted models","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"The package provides several methods that can be applied to fitted models. The names are similar to those in the GLM.jl package.","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"adjr2: adjusted R² for a linear model (an alias for adjr²)\naic: Akaike's Information Criterion\naicc: corrected Akaike's Information Criterion for small sample sizes\nbic: Bayesian Information Criterion\ncoef: estimates of the coefficients in the model\nconfint: confidence intervals for coefficients\ncooksdistance: Cook's distance for each observation\ndeviance: measure of the model fit, weighted residual sum of squares for lm's\ndispersion: dispersion (or scale) parameter for a model's distribution\ndof: number of degrees of freedom consumed in the model\ndof_residual: degrees of freedom for residuals, when meaningful\nfitted: fitted values of the model\nglm: fit a generalized linear model (an alias for fit(GeneralizedLinearModel, ...))\nlm: fit a linear model (an alias for fit(LinearModel, ...))\nloglikelihood: log-likelihood of the model\nmodelmatrix: design matrix\nnobs: number of rows, or sum of the weights when prior weights are specified\nnulldeviance: deviance of the model with all predictors removed\nnullloglikelihood: log-likelihood of the model with all predictors removed\npredict: predicted values of the dependent variable from the fitted model\nr2: R² of a linear model (an alias for r²)\nresiduals: vector of residuals from the fitted model\nresponse: model response (a.k.a the dependent variable)\nstderror: standard errors of the coefficients\nvcov: variance-covariance matrix of the coefficient estimates","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression","page":"Home","title":"FuzzyCRegression","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FuzzyCRegression.jl, which implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"theory.md\", \"package.md\", \"API.md\"]","category":"page"},{"location":"API/","page":"API","title":"API","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"API/#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.Model}","page":"API","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients for fcr model output\n\nArguments\n\nresults::Model Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\nplot:Bool Plot histogram of coefficients\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.fit-Tuple{}","page":"API","title":"FuzzyCRegression.fit","text":"fit()\n\nImplements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)\n\nArguments\n\ny::Vector dependent var\nx::Matrix variables with heterogeneous coefficients (defaults to vector of 1's for FEs)\nZ:matrix matrix of controls\nunit:Vector vector of unit IDs\nt::Vector time vector (optional)\nG::Integer number of groups (default = 2)\nm::Real fuzzy tuning parameter\nstartvals::Integer number of starting values (default = 1,000)\ncores::Integer number of cores (default = 1)\n\nReturns:\n\nstruct with the following methods:\n\ncoefficients::Vector: vector of coefficients\nweights::Matrix: group weights\nSE::Vector: standard errors (optional)\nobjective::Number value of objective function at minimum\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.information-Tuple{FuzzyCRegression.Model}","page":"API","title":"FuzzyCRegression.information","text":"information()\n\nCalculates inforation criteria for fcr model output, used for selecting optimal number of groups\n\nArguments\n\nresults::Model Model type from fcr output\ncriterion::String Information criterion (AIC, BIC, or HQC)\n\n\n\n\n\n","category":"method"}]
}
