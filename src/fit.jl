 # define FCRModel type for fitted model
struct FCRModel
    df
    y_name
    X_names
    Z_names
    coef
    value
    y
    X
    Z
    unit
    timed
    T
    N
    G
    m
end

 """
    fit()

Fits the FCR model

# Arguments
    - `df`: name of dataframe (if missing, data must be passed as arrays)
    - `y`: column name or array holding values of the dependent variable (required)
    - `X`: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (required)
    - `Z`: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients
    - `G`: number of groups (required)
    - `m`: regularization parameter (greater than 1), where group assignment becomes binary as m approaches 1 (default = 1.5)
    - `unit`: column name or array with unit identifier (if panel structure)
    - `time`: column name or array with time indicators (if panel structure)
    - `startvals`: number of starting values for the minimization routine (default = 100)
    - `cores`: number of parallel workers (default = 1)
"""
function fit(;df=nothing,y,X,Z=nothing,unit=nothing,t=nothing,G,m=1.5,startvals=10,parallel=false)
    if ≠(df,nothing)
        y_name = y
        X_names = X
        Z_names = Z
        y = Matrix(df[:,y])
        if ≠(X,nothing)
            X = Matrix(df[:,X])
        elseif ≠(Z,nothing)
            Z = Matrix(df[:,Z])
        elseif ≠(unit,nothing)
            unit = Matrix(df[:,unit])
        elseif ≠(t,nothing)
            t = Matrix(df[:,t])
        end
    else
        X_names = nothing
        y_name = nothing
        Z_names = nothing
    end

    if unit === nothing
        unit = collect(1:length(y))
    end
    #sort vectors
    y = y[sortperm(unit)]
    X = X[sortperm(unit),:]


    if Z === nothing 
        Z = Array{Any}(undef,length(y),0)
    else 
        Z = Z[sortperm(unit),:]
    end

    # case with no panel dimension
    if t === nothing 
        t = ones(length(y),1)
    else
        t = t[sortperm(unit)]
    end

    # matrix of time dummies
    T = length(unique(t))
    if T > 1
        timed = dummycreate(t)
    else 
        timed = ones(length(y),1)
    end

    N = ÷(length(y), T)

    function objective(a)
        #create matrix of error terms 
        m = max(m,1.1)
        amat = [reshape(a[1:G*T*(size(X,2))],(G,T*size(X,2))) repeat(a[G*T*(size(X,2))+1:end]',G)]'
        ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*amat,(T,N,G)),[2,1,3])
        ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))
        ϵ_g = repeat(ϵ,1,1,G)

        #weights 
        wgt = reshape(sum((reshape(ϵ_g[:,:,1],(N,1,G))./ϵ_g).^(1/(m-1)),dims=2).^(-1),(N,G))

        #value of objective function:
        return (sum(sum((wgt.^m).*ϵ,dims=2),dims=1)[1,1,1])/(N*T)
    end

    # homogeneous specification
    if T > 1
        baseline = ([ones(length(y)) Z timed]'*[ones(length(y)) Z timed])\([ones(length(y)) Z timed]'*y)
    else 
        baseline =([ones(length(y)) Z]'*[ones(length(y)) Z])\([ones(length(y)) Z]'*y)
    end

    # pre-allocate
    argmin = zeros(G*T*size(X,2)+size(Z,2),startvals)
    min = zeros(startvals,1)

    #loop over starting values
    for i=1:startvals
        if size(Z,2) > 0
            sval = [rand(G*T*size(X,2)); baseline[2:size(Z,2)+1]]
        else
            sval = rand(G*T*size(X,2))
        end
        results = optimize(objective, sval, LBFGS(),autodiff=:forward)   
        min[i] = results.minimum
        argmin[:,i] = results.minimizer
    end

    # select minimizing coefficients across starting values
    value = findmin(min)[1]
    coef = argmin[:,findmin(min)[2][1]]
    
    #return FCRModel type
    return FCRModel(df, y_name, X_names, Z_names,coef, value, y, X, Z, unit, timed, T, N, G, m)
end


