

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
    T = results.T
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

Standard erros from fitted model
"""
function stderror(results::FCRModel)
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N

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
    H = ForwardDiff.hessian(objective,coef)
    GRAD = fcrgradient(coef)
    # sandwich matrix
    V = GRAD*GRAD'
    # variance covariance matrix
    vcov = (inv(H./N)*(V./N)*inv(H./N))./N
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
    H = ForwardDiff.hessian(objective,coef)
    GRAD = fcrgradient(coef)
    # sandwich matrix
    V = GRAD*GRAD'
    # return variance covariance matrix
    return (inv(H./N)*(V./N)*inv(H./N))./N
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
    T = results.T
    timed = results.timed
    N = results.N
    coefs = results.coeff

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
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N
    coefs = results.coeff

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
    residuals()

Get the vector of residuals from the fitted model
"""
function residuals(results::FCRModel)
    y = results.y
    X = results.X
    Z = results.Z
    G = results.G
    T = results.T
    timed = results.timed
    N = results.N
    coefs = results.coeff

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
    distribution()

Calculates distribution of weighted coefficients from fitted model

# Arguments
- `results::Model` Model type from fcr output
- `index::Integer` Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)
"""
function distribution(results::FCRModel; index = 1)
    #grab data and coeffs
    G = results.G
    T = results.T
    coeffs = results.coef[G*T*(index-1)+1:G*T*(index)]

    #outer product of group coefficients with group weights
    return results.weight*coeffs
end


  