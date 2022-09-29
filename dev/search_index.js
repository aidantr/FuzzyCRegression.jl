var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FuzzyCRegression","category":"page"},{"location":"#FuzzyCRegression","page":"Home","title":"FuzzyCRegression","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FuzzyCRegression.jl, which implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FuzzyCRegression]","category":"page"},{"location":"#FuzzyCRegression.distribution-Tuple{FuzzyCRegression.Model}","page":"Home","title":"FuzzyCRegression.distribution","text":"distribution()\n\nCalculates distribution of weighted coefficients for fcr model output\n\nArguments\n\nresults::Model Model type from fcr output\nindex::Integer Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)\nplot:Bool Plot histogram of coefficients\n\n\n\n\n\n","category":"method"},{"location":"#FuzzyCRegression.fit-Tuple{}","page":"Home","title":"FuzzyCRegression.fit","text":"fit()\n\nImplements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)\n\nArguments\n\ny::Vector dependent var\nx::Matrix variables with heterogeneous coefficients (defaults to vector of 1's for FEs)\nZ:matrix matrix of controls\nunit:Vector vector of unit IDs\nt::Vector time vector (optional)\nG::Integer number of groups (default = 2)\nm::Real fuzzy tuning parameter\nstartvals::Integer number of starting values (default = 1,000)\ncores::Integer number of cores (default = 1)\n\nReturns:\n\nstruct with the following methods:\n\ncoefficients::Vector: vector of coefficients\nweights::Matrix: group weights\nSE::Vector: standard errors (optional)\nobjective::Number value of objective function at minimum\n\n\n\n\n\n","category":"method"},{"location":"#FuzzyCRegression.information-Tuple{FuzzyCRegression.Model}","page":"Home","title":"FuzzyCRegression.information","text":"information()\n\nCalculates inforation criteria for fcr model output, used for selecting optimal number of groups\n\nArguments\n\nresults::Model Model type from fcr output\ncriterion::String Information criterion (AIC, BIC, or HQC)\n\n\n\n\n\n","category":"method"}]
}
