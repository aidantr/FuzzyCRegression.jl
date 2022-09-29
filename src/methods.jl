

 """
    information()

Calculates inforation criteria for fcr model output, used for selecting optimal number of groups

# Arguments
- `results::Model` Model type from fcr output
- `criterion::String` Information criterion (AIC, BIC, or HQC)
"""
function information(results::Model; criterion = "AIC")
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
    if criterion == "AIC"
        return N*log(sum(sum(abs.(ϵ),dims=2),dims=1)[1,1]/N)+2*K

    elseif criterion == "BIC"
        return "need to fix BIC"

    elseif criterion == "HQC"
        return N*log(sum(sum(abs.(ϵ),dims=2),dims=1)[1,1]/N)+2*K*log(N)
    end
end


"""
    distribution()

Calculates distribution of weighted coefficients for fcr model output

# Arguments
- `results::Model` Model type from fcr output
- `index::Integer` Column index of variable in X matrix to calculate coefficient distibution for (defaults to 1)
- `plot:Bool` Plot histogram of coefficients
"""
function distribution(results::Model; index = 1)

    #grab data and coeffs
    G = results.G
    T = results.T
    coeffs = results.coef[G*T*(index-1)+1:G*T*(index)]

    #outer product of group coefficients with group weights
    return results.weight*coeffs

end


  