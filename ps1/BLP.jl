using Random, Distributions, LinearAlgebra

Random.seed!(2);

P = P_star
ns = 100

### (step -2) reshape data for GMM
Xc_1 = 0.5 .* ([X_2[3,:]; X_2[1,:]; X_2[2,:]] .+ [X_2[2,:]; X_2[3,:]; X_2[1,:]])
Xc_2 = 0.5 .* ([X_3[3,:]; X_3[1,:]; X_3[2,:]] .+ [X_3[2,:]; X_3[3,:]; X_3[1,:]])
Z_mat = [reshape(W, (J*M, 1)) reshape(Z, (J*M, 1)) Xc_1 Xc_2]
X = [reshape(P, (J*M,1)) reshape(X_1, (J*M,1)) reshape(X_2, (J*M,1)) reshape(X_3, (J*M,1))]

### (step -1) guess parameters
J   = 3
M   = 100
α_init   = 5
β₁_init  = 5
β₂_init  = 5
β₃_init  = 5
σ_α_init = 5


### (step 0) redraw individuals
ν_p  = rand(LogNormal(0,1), ns)


function BLP_1(α, β₁, β₂, β₃, σ_α, Φ_inv_old, t)
    ### (steps 1-2)
    A = α .+ σ_α*ν_p
    δ_star = ones(J,M)
    for m in 1:M


        function g(δm)
            δμ1m = δm[1] .- Α*P[1,m]
            δμ2m = δm[2] .- Α*P[1,m]
            δμ3m = δm[3] .- Α*P[1,m]
        
            den = 1 .+ exp.(δμ1m) .+ exp.(δμ2m) .+ exp.(δμ3m)
        
            sm = [mean((exp.(δμ1m)) ./ den),
                mean((exp.(δμ2m)) ./ den),
                mean((exp.(δμ3m)) ./ den)]
            sm = max.(1e-12, sm) # avoid taking log of 0

            return(δm + (log.(S[:,m]) - log.(sm)))
        end

        δ_star[:,m] = fixedpointmap(g; iv = ones(3), tolerance=1E-7)[1]
        
    end

    ### (step 3) error term and moment value
    ω = δ_star .- (-α .* P + X_1 .* β₁ + X_2 .* β₂ + X_3 .* β₃)
    ω = reshape(ω, (J*M,1))
    if (t == 1) # in step 1, use (Z'Z)^-1 for weighting matrix
        obj = transpose(ω) * Z_mat * inv(transpose(Z_mat) * Z_mat) * transpose(Z_mat) * ω
    else # otherwise use consistent weighting matrix
        obj = transpose(ω) * Z_mat * Φ_inv_old * transpose(Z_mat) * ω
    end

    ### (step 4) GMM Estimator
    if (t == 1) # in step 1, use (Z'Z)^-1 for weighting matrix
        θ_hat = inv(transpose(X) * Z_mat * inv(transpose(Z_mat) * Z_mat) * transpose(Z_mat) * X) * transpose(X) * Z_mat * inv(transpose(Z_mat) * Z_mat) * transpose(Z_mat) * reshape(δ_star, (J*M,1))
    else # in step t > 1, use weighting matrix calculated in step 1
        θ_hat = inv(transpose(X) * Z_mat * Φ_inv_old * transpose(Z_mat) * X) * transpose(X) * Z_mat * Φ_inv_old * transpose(Z_mat) * reshape(δ_star, (J*M,1))
    end

    α_hat  = θ_hat[1]
    β₁_hat = θ_hat[2]
    β₂_hat = θ_hat[3]
    β₃_hat = θ_hat[4]

    ω_hat =  δ_star .- (-α_hat .* P + X_1 .* β₁_hat + X_2 .* β₂_hat + X_3 .* β₃_hat)
    ω_hat = reshape(ω_hat,(J*M,1))
    if (t == 1) # in first iteration calculate weighting matrix
        Φ     = transpose(Z_mat) * ω_hat * transpose(ω_hat) * Z_mat
        Φ_inv = inv(Φ)
    else # otherwise, keep old weighting matrix
        Φ_inv = Φ_inv_old
    end
    t += 1
    return(α_hat, β₁_hat, β₂_hat, β₃_hat, Φ_inv, obj, t)
end

BLP_1(α_init, β₁_init, β₂_init, β₃_init, σ_α_init, I, 1)
