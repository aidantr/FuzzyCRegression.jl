"""
aic()

Calculates Aikike Inforation Criteria for fitted FCR model, used for selecting optimal number of groups
"""
function aic(results::FCRModel; criterion = "AIC")
    #grab data
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.Tx
    timed = results.timed
    N = results.N

    # matrix of error terms
    coefmat = [reshape(results.coef[1:G*T*(size(X,2))],(G,T*size(X,2))) repeat(results.coef[G*T*(size(X,2))+1:end]',G)]'
    ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*coefmat,(T,N,G)),[2,1,3])
    ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))

    # calculate criterion
    K = length(results.coef)
    return N*log(sum(sum(abs.(ϵ),dims=2),dims=1)[1,1]/N)+2*K
end

"""
bic()

Calculates Bayesian Inforation Criteria for fitted FCR model, used for selecting optimal number of groups
"""
function bic(results::FCRModel)
    #grab data
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N

    # matrix of error terms
    coefmat = [reshape(results.coef[1:G*T*(size(X,2))],(G,T*size(X,2))) repeat(results.coef[G*T*(size(X,2))+1:end]',G)]'
    ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*coefmat,(T,N,G)),[2,1,3])
    ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))

    # calculate criterion
    K = length(results.coef)
    return N*log(sum(sum(abs.(ϵ),dims=2),dims=1)[1,1]/N)+2*K*log(N)

end

"""
coef()

Extract coefficients from fitted model
"""
function coef(results::FCRModel)
    return results.coef
end

"""
stderror()

Standard errors from fitted model
"""
function stderror(results::FCRModel)
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N
    m = results.m

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

    function fcrgradient(a)
        function unit_gradient(a,n) 
            function unit_objective(a)
                #create matrix of error terms
                m = max(m,1.1) 
                amat = [reshape(a[1:G*T*(size(X,2))],(G,T*size(X,2)))  repeat(a[G*T*(size(X,2))+1:end]',G)]'
                ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*amat,(T,N,G)),[2,1,3])
                ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))
                ϵ_g = repeat(ϵ,1,1,G)

                #weights 
                wgt = reshape(sum((reshape(ϵ_g[:,:,1],(N,1,G))./ϵ_g).^(1/(m-1)),dims=2).^(-1),(N,G))

                #value of objective function for each unit
                (sum((wgt.^m).*ϵ,dims=2)[n,1,1])/(N*T)
            end
            return ForwardDiff.gradient(unit_objective,a)
        end
        GRAD = zeros(G*T*size(X,2)+size(Z,2),N)
        for i = 1:N
            GRAD[:,i] = unit_gradient(a,i)
        end
        return GRAD
    end

    # gradient and hessian at minimum
    H = ForwardDiff.hessian(objective,results.coef)
    GRAD = fcrgradient(results.coef)
    # sandwich matrix
    V = GRAD*GRAD'
    # variance covariance matrix
    vcov = (inv(H./N)*(V./N)*inv(H./N))
    # return standard errors
    return diag(vcov).^(1/2)
end

"""
vcov()

Variance covariance matrix from fitted model
"""
function vcov(results::FCRModel)
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N
    m = results.m

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

    function fcrgradient(a)
        function unit_gradient(a,n) 
            function unit_objective(a)
                #create matrix of error terms 
                m = max(m,1.1)
                amat = [reshape(a[1:G*T*(size(X,2))],(G,T*size(X,2)))  repeat(a[G*T*(size(X,2))+1:end]',G)]'
                ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*amat,(T,N,G)),[2,1,3])
                ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))
                ϵ_g = repeat(ϵ,1,1,G)

                #weights 
                wgt = reshape(sum((reshape(ϵ_g[:,:,1],(N,1,G))./ϵ_g).^(1/(m-1)),dims=2).^(-1),(N,G))

                #value of objective function for each unit
                (sum((wgt.^m).*ϵ,dims=2)[n,1,1])/(N*T)
            end
            return ForwardDiff.gradient(unit_objective,a)
        end
        GRAD = zeros(G*T*size(X,2)+size(Z,2),N)
        for i = 1:N
            GRAD[:,i] = unit_gradient(a,i)
        end
        return GRAD
    end

    # gradient and hessian at minimum
    H = ForwardDiff.hessian(objective,results.coef)
    GRAD = fcrgradient(results.coef)
    # sandwich matrix
    V = GRAD*GRAD'
    # return variance covariance matrix
    return (inv(H./N)*(V./N)*inv(H./N))
end

"""
weights()

Calculate group weights from fitted model, using modal group membership
"""
function weights(results::FCRModel)
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    m = results.m
    T = results.T
    timed = results.timed
    N = results.N
    coefs = results.coef

    #create matrix of error terms
    m = max(m,1.1) 
    coefmat = [reshape(coefs[1:G*T*(size(X,2))],(G,T*size(X,2)))  repeat(coefs[G*T*(size(X,2))+1:end]',G)]'
    ϵ = permutedims(reshape(repeat(y,1,G) - [X.*timed Z]*coefmat,(T,N,G)),[2,1,3])
    ϵ = reshape(sum(ϵ.^2,dims=2),(N,G))
    ϵ_g = repeat(ϵ,1,1,G)

    #weights 
    return reshape(sum((reshape(ϵ_g[:,:,1],(N,1,G))./ϵ_g).^(1/(m-1)),dims=2).^(-1),(N,G))
end

"""
predict()

Obtain predicted values of the dependent variable from the fitted model, using modal group membership
"""
function predict(results::FCRModel)
    return "In progress"
end

"""
residuals()

Get the vector of residuals from the fitted model
"""
function residuals(results::FCRModel)
    return "In progress"
end


"""
distribution()

Calculates distribution of weighted coefficients from fitted model

# Arguments
- `results::Model` Model type from fcr output
- `index::Integer` Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)
"""
function distribution(results::FCRModel; index = 1)
    G = results.G
    T = results.T
    coeffs = results.coef[G*T*(index-1)+1:G*T*(index)]
    wgts = weights(results)
    #outer product of group coefficients with group weights
    return wgts*coeffs
end

"""
    coefnames()

Returns names of coefficients from fitted model
"""
function coefnames(results::FCRModel)
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    X_names = results.X_names
    Z_names = results.Z_names
    coef_names = Array{String}(undef, size(X,2)*G*T+size(Z,2))
    for X_var = 1:size(X,2) 
        for g = 1:G
            if T>1
                for t = 1:T
                    if results.df === nothing
                        newname = string("X",X_var," (g=",g,", t=",t,")")
                    else
                        newname = string(X_names[X_var]," (g=",g,", t=",t,")")
                    end
                    coef_names[(X_var-1)*G*T+(g-1)*T+t] = newname
                end
            else
                if results.df === nothing
                    newname = string("X",X_var," (g=",g,")")
                else
                    newname = string(X_names[X_var]," (g=",g,")")
                end
                coef_names[(X_var-1)*G+g] = newname
            end
        end
    end
    if size(Z,2) > 0
        for Z_var = 1:size(Z,2)
            if results.df === nothing
                newname = string("Z",Z_var)
            else
                newname = string(Z_names[Z_var])
            end
            coef_names[size(X,2)*G*T+Z_var] = newname
        end
    end
    return coef_names
end



"""
    summarize()

Summarizes results from fitted model

"""
function summarize(results::FCRModel;level=0.95)
    cc = coef(results)
    se = stderror(results)
    tt = cc ./ se
    p = ccdf.(Ref(FDist(1, results.N-1)), abs2.(tt))
    ci = se*quantile(TDist(results.N-1), (1-level)/2)
    levstr = isinteger(level*100) ? string(Integer(level*100)) : string(level*100)

    ctf = TextFormat(
    up_right_corner = ' ',
    up_left_corner = ' ',
    bottom_left_corner=' ',
    bottom_right_corner= ' ',
    up_intersection= '─',
    left_intersection= ' ',
    right_intersection= ' ',
    middle_intersection= '─',
    bottom_intersection= '─',
    column= ' ',
    hlines=[ :begin, :header, :end])

    #return ct
    return pretty_table(hcat(coefnames(results),cc,se,tt,p,cc+ci,cc-ci); tf=ctf,header = ["","Estimate","Std. Error","t value","Pr(>|t|)","Lower $levstr%","Upper $levstr%"])
    
end 

"""
    confint()

Summarizes results from fitted model

"""
function confint(results::FCRModel;level=0.95)
    cc = coef(results)
    se = stderror(results)
    ci = se*quantile(TDist(results.N-1), (1-level)/2)
    lowerCI = cc-ci
    upperCI = cc+ci
    
    return [lowerCI upperCI]
end 
