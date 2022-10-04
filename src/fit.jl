
# define FCRModel type for fitted model
struct FCRModel
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

Implements fuzzy clustering regression estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)

# Arguments
- `y::Vector` dependent var
- `x::Matrix` variables with heterogeneous coefficients (defaults to vector of 1's for FEs)
- `Z:matrix` matrix of controls
- `unit:Vector` vector of unit IDs
- `t::Vector` time vector (optional)
- `G::Integer` number of groups (default = 2)
- `m::Real` fuzzy tuning parameter
- `startvals::Integer` number of starting values (default = 1,000)
- `cores::Integer` number of cores (default = 1)

# Returns:
 `struct` with the following methods:
- `coefficients::Vector`: vector of coefficients
- `weights::Matrix`: group weights
- `SE::Vector`: standard errors (optional)
- `objective::Number` value of objective function at minimum
"""
function fit(;y,unit=nothing, X=nothing, Z=nothing,t=nothing,G=2,m=1.1,startvals=10,parallel=false)
    #sort vectors
    y = y[sortperm(unit)]

    if X === nothing
        X = ones(length(y),1)
    else
        X = X[sortperm(unit),:]
    end

    if Z === nothing 
        Z = Array{Any}(undef,size(y),0)
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
        baseline = [ones(length(y)) Z timed]\y 
    else 
        baseline = [ones(length(y)) Z]\y 
    end

    # pre-allocate
    argmin = zeros(G*T*size(X,2)+size(Z,2),startvals)
    min = zeros(startvals,1)

    #loop over starting values
    for i=1:startvals
        sval = [rand(G*T*size(X,2)); baseline[2:size(Z,2)+1]]
        results = optimize(objective, sval, LBFGS(),autodiff=:forward)   
        min[i] = results.minimum
        argmin[:,i] = results.minimizer
    end

    # select minimizing coefficients across starting values
    value = findmin(min)[1]
    coef = argmin[:,findmin(min)[2][1]]
    
    #return Model struct
    return FCRModel(coef, value, y, X, Z, unit, timed, T, N, G, m)
end


