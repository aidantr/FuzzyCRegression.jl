var documenterSearchIndex = {"docs":
[{"location":"theory/#Fuzzy-C-Regression","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of heterogeneity. It extends the \"Fuzzy C-Means\" clustering algorithm to regression settings.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"FCR can be used to estimate \"grouped fixed effects\" (a constant term for each group) as well as heterogeneous coefficients.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"In panel settings, each unit belongs to the same group over time, meaning that separate coefficients are estimated for each group-time pair.","category":"page"},{"location":"theory/#FCR-objective-function","page":"Fuzzy C-Regression","title":"FCR objective function","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"y = sum_g=1^Gmu_g theta_gX+varepsilon","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where X are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g represent group weights. FCR is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"The FCR objective function takes the form:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g xrightVert ^2right","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"mu_gleft(yxthetamright)=left(sum_h=1^GfracleftVert y-theta_g xrightVert ^2left(m-1right)leftVert y-theta_h xrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"J^FCRleft(thetaright)=Eleftleft(sum_g=1^GleftVert y-theta_g xrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"theory/#Useful-properties","page":"Fuzzy C-Regression","title":"Useful properties","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Fast: FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.\nCustomizable: The “fuzziness” of the FCR groups is governed by the regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.\nInference: Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors. ","category":"page"},{"location":"package/#FuzzyCRegression.jl-tutorial","page":"Examples","title":"FuzzyCRegression.jl tutorial","text":"","category":"section"},{"location":"package/#Installation","page":"Examples","title":"Installation","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"using Pkg; Pkg.add(\"FuzzyCRegression\")","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"will install the package and its dependencies, which include Optim.jl for minimization and ForwardDiff.jl for automatic differentiation.","category":"page"},{"location":"package/#Fitting-the-FCR-model","page":"Examples","title":"Fitting the FCR model","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"There are two ways to fit an FCR model, depending on whether the data is stored as a DataFrame or as a set of arrays.","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"If the dataset is stored as a DataFrame, the model can be fit using fit(df,y,X,G,m,...), where the variables are referenced by their column names. G and m specify the number of groups and regularization parameter, respectively. For example, using the iris dataset from RDatasets:","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"using FuzzyCRegression, RDatasets\n\niris = dataset(\"datasets\", \"iris\")\n\nfcr_model = fit(df = iris,y = \"SepalLength\", x = [\"SepalWidth\", \"PetalWidth\"], G = 3, m = 1.5)\nsummarize(fcr_model)","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"An advantage of this approach is that the estimated coefficients are labeled by variable name.","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"Alternatively, the data can be passed directly as arrays:","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"y = iris.SepalLength\nX = [iris.SepalWidth iris.PetalWidth]\n\nfcr_model = fit(y, X, G = 3, m = 1.5)\nsummarize(fcr_model)","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"The arguments for fitting the model are:","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"df: name of dataframe (if missing, data must be passed as arrays)\ny: column name or array holding values of the dependent variable\nX: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (defaults to a constant)\nZ: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients\nG: number of groups\nm: regularization parameter (greater than 1), where group assignment becomes binary as m rightarrow 1\nunit: column name or array with unit identifier (if panel structure)\ntime: column name or array with time indicators (if panel structure)\nstartvals: number of starting values for the minimization routine (default = 10)\ncores: number of parallel workers (default = 1)","category":"page"},{"location":"package/#Methods-applied-to-fitted-models","page":"Examples","title":"Methods applied to fitted models","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"The package provides several methods that can be applied to fitted models. The names are similar to those in the GLM.jl package. Full documentation for these methods can be found here.","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"aic: Akaike's Information Criterion\nbic: Bayesian Information Criterion\ncoef: estimates of the coefficients in the model\nconfint: confidence intervals for coefficients\ndistribution: distribution of coefficients using group weights\npredict: obtain predicted values of the dependent variable from the fitted model, using modal group membership\nresiduals: vector of residuals from the fitted model, using modal group membership\nstderror: standard errors of the coefficients\nsummarize: summarize model results\nvcov: variance-covariance matrix of the coefficient estimates","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"For example, to plot the distribution of coefficients on SepalWidth from the fitted model above:","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"SepalWidth_coefs = distribution(fcr_model,\"SepalWidth\")\n\nusing Gadfly\nplot(SepalWidth_coefs, Geom.hist)","category":"page"},{"location":"package/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"package/#Simple-example","page":"Examples","title":"Simple example","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"To illustrate the package's functionality, we start with a simple example estimating grouped fixed effects with no controls. ","category":"page"},{"location":"package/#More-complicated-example","page":"Examples","title":"More complicated example","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"We now turn to a more complicated example where the data has a panel structure, we estimate heterogeneous coefficients on the independent variable for each time period, and we include a common control variable. ","category":"page"},{"location":"package/#Choosing-G","page":"Examples","title":"Choosing G","text":"","category":"section"},{"location":"package/#Choosing-m","page":"Examples","title":"Choosing m","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression.jl-manual","page":"Home","title":"FuzzyCRegression.jl manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022) in Julia. ","category":"page"},{"location":"#Fuzzy-C-Regression-estimator","page":"Home","title":"Fuzzy C-Regression estimator","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of unobserved heterogeneity. It extends the \"Fuzzy C-Means\" clustering algorithm to regression settings.","category":"page"},{"location":"","page":"Home","title":"Home","text":"FCR can be used to estimate \"grouped fixed effects\" (a constant term for each group) as well as heterogeneous coefficients.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In panel settings, each unit belongs to the same group over time, meaning that separate coefficients are estimated for each group-time pair.","category":"page"},{"location":"#FCR-objective-function","page":"Home","title":"FCR objective function","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"","page":"Home","title":"Home","text":"y = sum_g=1^Gmu_g theta_gX+varepsilon","category":"page"},{"location":"","page":"Home","title":"Home","text":"where X are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g represent group weights which sum to 1. FCR is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The FCR objective function takes the form:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g XrightVert ^2right","category":"page"},{"location":"","page":"Home","title":"Home","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"","page":"Home","title":"Home","text":"mu_gleft(yXthetamright)=left(sum_h=1^GfracleftVert y-theta_g XrightVert ^2left(m-1right)leftVert y-theta_h XrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"","page":"Home","title":"Home","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"","page":"Home","title":"Home","text":"J^FCRleft(thetaright)=Eleftleft(sum_g=1^GleftVert y-theta_g XrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"","page":"Home","title":"Home","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"#Useful-properties","page":"Home","title":"Useful properties","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fast: FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.\nCustomizable: The “fuzziness” of the FCR groups is governed by the regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.\nInference: Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors.","category":"page"},{"location":"#Package-installation","page":"Home","title":"Package installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FuzzyCRegressions.jl can be installed using:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pkg.add(\"FuzzyCRegression\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"This adds the latest version of the package and its dependencies, which include Optim.jl for minimization and ForwardDiff.jl for automatic differentiation. ","category":"page"},{"location":"#Fitting-the-FCR-model","page":"Home","title":"Fitting the FCR model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are two ways to fit an FCR model, depending on whether the data is stored as a DataFrame or as a set of arrays.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If the dataset is stored as a DataFrame, the model can be fit using fit(df,y,X,G,m,...), where the variables are referenced by their column names. G specifies the number of groups and m sets the regularization parameter, where group membership becomes binary as m rightarrow 1^+ (tips for selecting these options are discussed below). For example, using the iris dataset from RDatasets:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FuzzyCRegression, RDatasets\n\niris = dataset(\"datasets\", \"iris\")\n\nfcr_model = fit(df = iris,y = \"SepalLength\", x = [\"SepalWidth\", \"PetalWidth\"], G = 3, m = 1.5)\nsummarize(fcr_model)","category":"page"},{"location":"","page":"Home","title":"Home","text":"An advantage of this approach is that the estimated coefficients are labeled by variable name.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Alternatively, the data can be passed directly as arrays:","category":"page"},{"location":"","page":"Home","title":"Home","text":"y = iris.SepalLength\nX = [iris.SepalWidth iris.PetalWidth]\n\nfcr_model = fit(y, X, G = 3, m = 1.5)\nsummarize(fcr_model)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The arguments for fitting the model are:","category":"page"},{"location":"","page":"Home","title":"Home","text":"df: name of dataframe (if missing, data must be passed as arrays)\ny: column name or array holding values of the dependent variable\nX: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (defaults to a constant)\nZ: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients\nG: number of groups\nm: regularization parameter (greater than 1), where group assignment becomes binary as m rightarrow 1\nunit: column name or array with unit identifier (if panel structure)\ntime: column name or array with time indicators (if panel structure)\nstartvals: number of starting values for the minimization routine (default = 10)\ncores: number of parallel workers (default = 1)","category":"page"},{"location":"#Methods-applied-to-fitted-models","page":"Home","title":"Methods applied to fitted models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package provides several methods that can be applied to fitted models. The names are similar to those in GLM.jl. Full documentation for these methods can be found here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"aic: Akaike's Information Criterion\nbic: Bayesian Information Criterion\ncoef: estimates of the coefficients in the model\nconfint: confidence intervals for coefficients\ndistribution: distribution of coefficients using group weights\npredict: obtain predicted values of the dependent variable from the fitted model, using modal group membership\nresiduals: vector of residuals from the fitted model, using modal group membership\nstderror: standard errors of the coefficients\nsummarize: summarize model results\nvcov: variance-covariance matrix of the coefficient estimates","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, to plot the distribution of coefficients on SepalWidth from the fitted model above:","category":"page"},{"location":"","page":"Home","title":"Home","text":"SepalWidth_coefs = distribution(fcr_model,\"SepalWidth\")\n\nusing Gadfly\nplot(SepalWidth_coefs, Geom.hist)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Choosing-G","page":"Home","title":"Choosing G","text":"","category":"section"},{"location":"#Choosing-m","page":"Home","title":"Choosing m","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"API/#FuzzyCRegression.aic-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.aic","text":"aic()\n\nCalculates Aikike Inforation Criteria for fitted FCR model, used for selecting optimal number of groups\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.bic-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.bic","text":"bic()\n\nCalculates Bayesian Inforation Criteria for fitted FCR model, used for selecting optimal number of groups\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.coef-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.coef","text":"coef()\n\nExtract coefficients from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients from fitted model\n\nArguments\n\nresults::Model Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.fit-Tuple{}","page":"API","title":"FuzzyCRegression.fit","text":"fit()\n\nImplements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)\n\nArguments\n\ny::Vector dependent var\nx::Matrix variables with heterogeneous coefficients (defaults to vector of 1's for FEs)\nZ:matrix matrix of controls\nunit:Vector vector of unit IDs\nt::Vector time vector (optional)\nG::Integer number of groups (default = 2)\nm::Real fuzzy tuning parameter\nstartvals::Integer number of starting values (default = 1,000)\ncores::Integer number of cores (default = 1)\n\nReturns:\n\nstruct with the following methods:\n\ncoefficients::Vector: vector of coefficients\nweights::Matrix: group weights\nSE::Vector: standard errors (optional)\nobjective::Number value of objective function at minimum\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.predict-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.predict","text":"predict()\n\nObtain predicted values of the dependent variable from the fitted model, using modal group membership\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.residuals-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.residuals","text":"residuals()\n\nGet the vector of residuals from the fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.stderror-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.stderror","text":"stderror()\n\nStandard errors from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.vcov-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.vcov","text":"vcov()\n\nVariance covariance matrix from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.weights-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.weights","text":"weights()\n\nCalculate group weights from fitted model, using modal group membership\n\n\n\n\n\n","category":"method"}]
}
