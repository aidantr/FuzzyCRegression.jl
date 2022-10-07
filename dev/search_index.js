var documenterSearchIndex = {"docs":
[{"location":"theory/#Fuzzy-C-Regression","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of heterogeneity. It extends the \"Fuzzy C-Means\" clustering algorithm to regression settings.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"FCR can be used to estimate \"grouped fixed effects\" (a constant term for each group) as well as heterogeneous coefficients.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"In panel settings, each unit belongs to the same group over time, meaning that separate coefficients are estimated for each group-time pair.","category":"page"},{"location":"theory/#FCR-objective-function","page":"Fuzzy C-Regression","title":"FCR objective function","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"y = sum_g=1^Gmu_g theta_gX+varepsilon","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where X are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g represent group weights. FCR is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"The FCR objective function takes the form:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"L^FCR_mleft(thetamuright)=Eleftsum_g=1^Gmu_g^mleftVert y-theta_g xrightVert ^2right","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"where m  1 is the regularization parameter. The weights are defined as","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"mu_gleft(yxthetamright)=left(sum_h=1^GfracleftVert y-theta_g xrightVert ^2left(m-1right)leftVert y-theta_h xrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"J^FCRleft(thetaright)=Eleftleft(sum_g=1^GleftVert y-theta_g xrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"theory/#Useful-properties","page":"Fuzzy C-Regression","title":"Useful properties","text":"","category":"section"},{"location":"theory/","page":"Fuzzy C-Regression","title":"Fuzzy C-Regression","text":"Fast: FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.\nCustomizable: The “fuzziness” of the FCR groups is governed by the regularization parameter m, where group assignment becomes binary as m rightarrow 1. Choosing m allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.\nInference: Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors. ","category":"page"},{"location":"package/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"package/","page":"Examples","title":"Examples","text":"We now illustrate the FCR estimator through two examples.","category":"page"},{"location":"package/#Example-1:","page":"Examples","title":"Example 1:","text":"","category":"section"},{"location":"package/#Example-2:","page":"Examples","title":"Example 2:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression.jl-Manual","page":"Home","title":"FuzzyCRegression.jl Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022) in Julia. The source code can be found here.","category":"page"},{"location":"#Fuzzy-C-Regression","page":"Home","title":"Fuzzy C-Regression","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of unobserved heterogeneity. It extends the \"Fuzzy C-Means\" clustering algorithm to regression.","category":"page"},{"location":"","page":"Home","title":"Home","text":"FCR can be used to estimate \"grouped fixed effects\" (a constant term for each group) as well as heterogeneous coefficients.","category":"page"},{"location":"#FCR-objective-function","page":"Home","title":"FCR objective function","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Consider a linear model with grouped heterogeneity:","category":"page"},{"location":"","page":"Home","title":"Home","text":"y = sum_g=1^Gmu_g theta_gX+varepsilon","category":"page"},{"location":"","page":"Home","title":"Home","text":"where X are covariates (which could simply be a constant), theta_g represent group-specific coefficients for groups g=1ldotsG, and mu_g represent group weights which sum to 1. FCR is concerned with jointly estimating theta and mu for each group.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The FCR objective function takes the form:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L^FCR_mleft(thetamuright)=mathbbEleftsum_g=1^Gmu_g^mleftVert y-theta_g XrightVert ^2right","category":"page"},{"location":"","page":"Home","title":"Home","text":"The regularization parameter m  1 governs the \"fuzziness\" of the FCR clusters, where group membership becomes binary as m rightarrow 1^+. The weights are defined as","category":"page"},{"location":"","page":"Home","title":"Home","text":"mu_gleft(yXthetamright)=left(sum_h=1^GfracleftVert y-theta_g XrightVert ^2left(m-1right)leftVert y-theta_h XrightVert ^2left(m-1right)right)^-1 text for  g=1ldotsG","category":"page"},{"location":"","page":"Home","title":"Home","text":"Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors y-theta_gX:","category":"page"},{"location":"","page":"Home","title":"Home","text":"J^FCRleft(thetaright)=mathbbEleftleft(sum_g=1^GleftVert y-theta_g XrightVert ^-2left(m-1right)right)^1-mright","category":"page"},{"location":"","page":"Home","title":"Home","text":"Thus, for fixed m, the FCR function is differentiable and can be written as a standard GMM problem. ","category":"page"},{"location":"#Useful-properties","page":"Home","title":"Useful properties","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fast: FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units.\nCustomizable: Choosing the fuzzy tuning parameter m allows FCR to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects when heterogeneity is not fully discrete.\nInference: Since FCR is a GMM problem, analytic standard errors can be derived from standard theory.","category":"page"},{"location":"#Package-installation","page":"Home","title":"Package installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FuzzyCRegressions.jl can be installed using:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pkg.add(\"FuzzyCRegression\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"This adds the latest version of the package and its dependencies, which include Optim.jl for minimization and ForwardDiff.jl for automatic differentiation. ","category":"page"},{"location":"#Fitting-the-FCR-model","page":"Home","title":"Fitting the FCR model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are two ways to fit an FCR model, depending on whether the data is stored as a DataFrame or as a set of arrays.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If the dataset is stored as a DataFrame, the model can be fit using fit(df, y, X, G, m, ...), where the variables are referenced by their column names. G specifies the number of groups and m sets the regularization parameter (tips for selecting these options are discussed below). fit returns an object of type FCRModel, to which a number of methods can be applied.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, using the iris dataset from RDatasets:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FuzzyCRegression, RDatasets\n\niris = dataset(\"datasets\", \"iris\")\n\nfcr_model = fit(df=iris, y=[\"PetalLength\"], X=[\"SepalWidth\",\"PetalWidth\"], G=3, m=1.5)\nsummarize(fcr_model)\n\n ────────────────────────────────────────────────────────────────────────────────────────────\n                      Estimate   Std. Error    t value     Pr(>|t|)   Lower 95%    Upper 95% \n ────────────────────────────────────────────────────────────────────────────────────────────\n   SepalWidth (g=1)   -1.85059     0.914888   -2.02275    0.0448874    -3.65842   -0.0427603\n   SepalWidth (g=2)   -1.89972     0.581725   -3.26568   0.00135507    -3.04922    -0.750229\n   SepalWidth (g=3)   0.867675       2.8992   0.299281     0.765143    -4.86118      6.59653\n  SepalLength (g=1)     1.5167     0.397579    3.81484   0.00019913     0.73108      2.30232\n  SepalLength (g=2)    1.70667     0.338777    5.03774   1.34499e-6     1.03724       2.3761\n  SepalLength (g=3)   0.412295       1.4735   0.279807     0.780014    -2.49936      3.32395\n ────────────────────────────────────────────────────────────────────────────────────────────","category":"page"},{"location":"","page":"Home","title":"Home","text":"An advantage of this approach is that the estimated coefficients in the regression output are labeled by variable name.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Alternatively, the data can be passed directly as arrays:","category":"page"},{"location":"","page":"Home","title":"Home","text":"fcr_model = fit(y=iris.PetalLength, X=[iris.SepalWidth iris.SepalLength], G=3, m=1.5)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The arguments for fitting the model are:","category":"page"},{"location":"","page":"Home","title":"Home","text":"df: name of dataframe (if missing, data must be passed as arrays)\ny: column name or array holding values of the dependent variable (required)\nX: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (required)\nZ: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients\nG: number of groups (required)\nm: regularization parameter (greater than 1), where group assignment becomes binary as m rightarrow 1 (default = 1.5)\nunit: column name or array with unit identifier (if panel structure)\ntime: column name or array with time indicators (if panel structure)\nstartvals: number of starting values for the minimization routine (default = 10)\ncores: number of parallel workers (default = 1)","category":"page"},{"location":"#Methods-applied-to-fitted-models","page":"Home","title":"Methods applied to fitted models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package provides several methods that can be applied to fitted models. The names are similar to those in GLM.jl and R. Full documentation for these functions can be found here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"aic: Akaike's Information Criterion\nbic: Bayesian Information Criterion\ncoef: estimates of the coefficients in the model\nconfint: confidence intervals for coefficients\ndistribution: distribution of coefficients using group weights\npredict: obtain predicted values of the dependent variable from the fitted model, using modal group membership\nresiduals: vector of residuals from the fitted model, using modal group membership\nstderror: standard errors of the coefficients\nsummarize: summarize model results\nvcov: variance-covariance matrix of the coefficient estimates","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, to plot the group-weighted distribution of coefficients on SepalWidth from the fitted model above:","category":"page"},{"location":"","page":"Home","title":"Home","text":"SepalWidth_coefs = distribution(fcr_model,\"SepalWidth\")\n\nusing Gadfly\nplot(SepalWidth_coefs, Geom.hist)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Choosing-G","page":"Home","title":"Choosing G","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package provides several data-driven approaches for choosing the number of groups. In particular, the aicand bic methods calculate information criteria that trade off the fit of the model against the number of parameters. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, continuing with the iris dataset, the minimum of the AIC and BIC indicate that the optimal group number is 3 or 4 according to these criteria.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FuzzyCRegression, Gadfly\n\nIC = zeros(15,2)\nfor g = 1:15\n    fcr_model = fit(df=iris, y=[\"SepalLength\"], X=[\"SepalWidth\",\"PetalWidth\"], G=g, m=1.5)\n    IC[g,1] = aic(fcr_model)\n    IC[g,2] = bic(fcr_model)\nend\nIC_norm = IC./sqrt(sum(IC.^2,1))\n\nplot(x=collect(1:15), y=IC_norm, Geom.point, Geom.line)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Choosing-m","page":"Home","title":"Choosing m","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The regularization parameter m determines the \"fuzziness\" of the FCR groups. Fuzzy clusters are helpful for two reasons:","category":"page"},{"location":"","page":"Home","title":"Home","text":"First, even in settings with discrete unobserved heterogeneity, noise means that group membership cannot be ascertained with certainty. Thus,             probabilistic clustering improves performance.\nSecond, in many applications heterogeneity may be continuous, and the fuzzy clusters allow to FCR to approximate its distribution.","category":"page"},{"location":"API/","page":"API","title":"API","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"API/#FuzzyCRegression.aic-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.aic","text":"aic()\n\nCalculates Aikike Inforation Criteria for fitted FCR model, used for selecting optimal number of groups\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.bic-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.bic","text":"bic()\n\nCalculates Bayesian Inforation Criteria for fitted FCR model, used for selecting optimal number of groups\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.coef-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.coef","text":"coef()\n\nExtract coefficients from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.coefnames-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.coefnames","text":"coefnames()\n\nReturns names of coefficients from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.confint-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.confint","text":"confint()\n\nReturns lower and upper confidence interval for fitted model, for specified significance level (default = 0.95)\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients from fitted model\n\nArguments\n\nresults::FCRModel Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.fit-Tuple{}","page":"API","title":"FuzzyCRegression.fit","text":"fit()\n\nFits the FCR model\n\nArguments\n\n- `df`: name of dataframe (if missing, data must be passed as arrays)\n- `y`: column name or array holding values of the dependent variable (required)\n- `X`: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (required)\n- `Z`: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients\n- `G`: number of groups (required)\n- `m`: regularization parameter (greater than 1), where group assignment becomes binary as m approaches 1 (default = 1.5)\n- `unit`: column name or array with unit identifier (if panel structure)\n- `time`: column name or array with time indicators (if panel structure)\n- `startvals`: number of starting values for the minimization routine (default = 100)\n- `cores`: number of parallel workers (default = 1)\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.predict-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.predict","text":"predict()\n\nObtain predicted values of the dependent variable from the fitted model, using modal group membership\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.residuals-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.residuals","text":"residuals()\n\nGet the vector of residuals from the fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.stderror-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.stderror","text":"stderror()\n\nStandard errors from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.summarize-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.summarize","text":"summarize()\n\nSummarizes results from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.vcov-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.vcov","text":"vcov()\n\nVariance covariance matrix from fitted model\n\n\n\n\n\n","category":"method"},{"location":"API/#FuzzyCRegression.weights-Tuple{FuzzyCRegression.FCRModel}","page":"API","title":"FuzzyCRegression.weights","text":"weights()\n\nCalculate group weights from fitted model, using modal group membership\n\n\n\n\n\n","category":"method"}]
}
