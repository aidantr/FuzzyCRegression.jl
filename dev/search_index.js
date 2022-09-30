var documenterSearchIndex = {"docs":
[{"location":"theory/#Fuzzy-Clustering-Regression","page":"FCR","title":"Fuzzy Clustering Regression","text":"","category":"section"},{"location":"theory/#FCR-Objective","page":"FCR","title":"FCR Objective","text":"","category":"section"},{"location":"theory/","page":"FCR","title":"FCR","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"y_it = sum_g=1^Gmu_g(i) theta_gX_i+varepsilon_it","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"where X_it are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g(i) represent group weights for unit i. Fuzzy Clustering Regression (FCR) is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"The FCR objective function takes the form:","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g xrightVert ^2right","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"mu_gleft(yxthetamright)=left(sum_h=1^GfracleftVert y-theta_g xrightVert ^2left(m-1right)leftVert y-theta_h xrightVert ^2left(m-1right)right)^-1g=1ldotsG","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y=theta_gX. Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. This has two important implications:","category":"page"},{"location":"theory/#Useful-Properties","page":"FCR","title":"Useful Properties","text":"","category":"section"},{"location":"theory/","page":"FCR","title":"FCR","text":"FRC can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units into groups. While FCR remains quite efficient with large datasets, computation time for iterative approaches becomes prohibitive.\nThe “fuzziness” of the FCR groups is governed by the","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows one to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it FCR can then recover the full distribution of effects.","category":"page"},{"location":"theory/","page":"FCR","title":"FCR","text":"Since FCR is a GMM problem, it's asymptotic properties follow from standard arguments and we are able to derive analytic standard errors. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression","page":"Home","title":"FuzzyCRegression","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FuzzyCRegression.jl, which implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"theory.md\", \"sparse_arrays_jacs.md\", \"numerical_algorithms.md\", \"example.md\", \"caching.md\", \"diagnostics.md\", \"tips.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.Model}","page":"Home","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients for fcr model output\n\nArguments\n\nresults::Model Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\nplot:Bool Plot histogram of coefficients\n\n\n\n\n\n","category":"method"},{"location":"#FuzzyCRegression.fit-Tuple{}","page":"Home","title":"FuzzyCRegression.fit","text":"fit()\n\nImplements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)\n\nArguments\n\ny::Vector dependent var\nx::Matrix variables with heterogeneous coefficients (defaults to vector of 1's for FEs)\nZ:matrix matrix of controls\nunit:Vector vector of unit IDs\nt::Vector time vector (optional)\nG::Integer number of groups (default = 2)\nm::Real fuzzy tuning parameter\nstartvals::Integer number of starting values (default = 1,000)\ncores::Integer number of cores (default = 1)\n\nReturns:\n\nstruct with the following methods:\n\ncoefficients::Vector: vector of coefficients\nweights::Matrix: group weights\nSE::Vector: standard errors (optional)\nobjective::Number value of objective function at minimum\n\n\n\n\n\n","category":"method"},{"location":"#FuzzyCRegression.information-Tuple{FuzzyCRegression.Model}","page":"Home","title":"FuzzyCRegression.information","text":"information()\n\nCalculates inforation criteria for fcr model output, used for selecting optimal number of groups\n\nArguments\n\nresults::Model Model type from fcr output\ncriterion::String Information criterion (AIC, BIC, or HQC)\n\n\n\n\n\n","category":"method"}]
}
