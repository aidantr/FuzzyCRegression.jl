var documenterSearchIndex = {"docs":
[{"location":"theory/#Fuzzy-C-Regression","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Fuzzy C-Regression, or FCR, is a method for estimating heterogeneous effects in settings with grouped patterns of heterogeneity. It extends the \"Fuzzy C-Means\" clustering algorithm to regression settings.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"FCR can be used to estimate \"grouped fixed effects\" (a constant term for each group) as well as heterogeneous coefficients on independent variables. ","category":"page"},{"location":"theory/#FCR-objective-function","page":"Fuzzy C-Regression","title":"FCR objective function","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"y_it = sum_g=1^Gmu_g(i) theta_gtX_it+varepsilon_it","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where X_it are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g(i) represent group weights for unit i. FCR is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"The FCR objective function takes the form:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g xrightVert ^2right","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"mu_gleft(yxthetamright)=left(sum_h=1^GfracleftVert y-theta_g xrightVert ^2left(m-1right)leftVert y-theta_h xrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"J^FCRleft(thetaright)=Eleftleft(sum_g=1^GleftVert y-theta_g xrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"theory/#Useful-properties","page":"Fuzzy C-Regression","title":"Useful properties","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"___Fast:___ FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.\n___Customizable:___ The “fuzziness” of the FCR groups is governed by the regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.\n___Inference:___ Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors. ","category":"page"},{"location":"package/#FuzzyCRegression-Tutorial","page":"Julia Implementation","title":"FuzzyCRegression Tutorial","text":"","category":"section"},{"location":"package/#Installation","page":"Julia Implementation","title":"Installation","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"Pkg.add(\"FuzzyCRegression\")","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"will install the package and its dependencies, which include Optim.jl for minimization and ForwardDiff.jl for automatic differentiation.","category":"page"},{"location":"package/#Fitting-the-FCR-model","page":"Julia Implementation","title":"Fitting the FCR model","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"There are two ways to fit an FCR model, using DataFrames or arrays. ","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"If the data is stored as a DataFrame, the model can be fit using fit(df,y,X,G,m,...), where the variables are referenced by their column names. For example, using the iris dataset from RDatasets:","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"using FuzzyCRegression, RDatasets\n\niris = dataset(\"datasets\", \"iris\")\n\nfcr_model = fit(iris,y = \"SepalWidth\", x = [\"1\",\"SepalLength\"], G=3, m=1.5)\nsummarize(fcr_model)","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"where \"1\" specifies a constant term. An advantage of this approach is that the estimated coefficients are labeled by variable name.","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"Alternatively, the data can be passed directly as arrays:","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"y = iris.SepalWidth\nX = [ones(length(y)) iris.SepalLength]\n\nfcr_model = fit(y, X, G=3, m=1.5)\nsummarize(fcr_model)","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"The optional arguments for fitting the model are:","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"Z: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients\nG: number of groups\nm: regularization parameter (greater than 1), where group assignment becomes binary as m rightarrow 1\nunit: unit identifier (if panel structure)\ntime: vector of time indicators (if panel structure)\nstartvals: number of starting values for the minimization routine (default = 100)\ncores: number of parallel workers (default = 1)","category":"page"},{"location":"package/#Methods-applied-to-fitted-models","page":"Julia Implementation","title":"Methods applied to fitted models","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"The package provides several methods that can be applied to fitted models. The names are similar to those in the GLM.jl package.","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"aic: Akaike's Information Criterion\nbic: Bayesian Information Criterion\ncoef: estimates of the coefficients in the model\nconfint: confidence intervals for coefficients\ndistribution: distribution of coefficients using group weights\npredict: obtain predicted values of the dependent variable from the fitted model, using modal group membership\nresiduals: vector of residuals from the fitted model, using modal group membership\nstderror: standard errors of the coefficients\nsummarize: summarize model results\nvcov: variance-covariance matrix of the coefficient estimates","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"For example, to plot the distribution of coefficients on SepalLength from the previous example:","category":"page"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"using Gadfly\n\nSepalLength_coefs = distribution(fcr_model,\"SepalLength\")\nplot(SepalLength_coefs, geom.hist)\n","category":"page"},{"location":"package/#Simple-example","page":"Julia Implementation","title":"Simple example","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"To illustrate the package's functionality, we start with a simple example estimating grouped fixed effects with no controls. ","category":"page"},{"location":"package/#More-complicated-example","page":"Julia Implementation","title":"More complicated example","text":"","category":"section"},{"location":"package/","page":"Julia Implementation","title":"Julia Implementation","text":"We now turn to a more complicated example where the data has a panel structure, we estimate heterogeneous coefficients on the independent variable for each time period, and we include a common control variable. ","category":"page"},{"location":"package/#Choosing-G","page":"Julia Implementation","title":"Choosing G","text":"","category":"section"},{"location":"package/#Choosing-m","page":"Julia Implementation","title":"Choosing m","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression.jl","page":"Home","title":"FuzzyCRegression.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FuzzyCRegression.jl, which implements the heterogeneous effects estimator of Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"theory.md\", \"package.md\", \"API.md\"]","category":"page"},{"location":"API/","page":"API","title":"API","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"API/#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.Model}","page":"API","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients for fcr model output\n\nArguments\n\nresults::Model Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\nplot:Bool Plot histogram of coefficients\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.fit-Tuple{}","page":"API","title":"FuzzyCRegression.fit","text":"fit()\n\nImplements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)\n\nArguments\n\ny::Vector dependent var\nx::Matrix variables with heterogeneous coefficients (defaults to vector of 1's for FEs)\nZ:matrix matrix of controls\nunit:Vector vector of unit IDs\nt::Vector time vector (optional)\nG::Integer number of groups (default = 2)\nm::Real fuzzy tuning parameter\nstartvals::Integer number of starting values (default = 1,000)\ncores::Integer number of cores (default = 1)\n\nReturns:\n\nstruct with the following methods:\n\ncoefficients::Vector: vector of coefficients\nweights::Matrix: group weights\nSE::Vector: standard errors (optional)\nobjective::Number value of objective function at minimum\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.information-Tuple{FuzzyCRegression.Model}","page":"API","title":"FuzzyCRegression.information","text":"information()\n\nCalculates inforation criteria for fcr model output, used for selecting optimal number of groups\n\nArguments\n\nresults::Model Model type from fcr output\ncriterion::String Information criterion (AIC, BIC, or HQC)\n\n\n\n\n\n","category":"method"}]
}
