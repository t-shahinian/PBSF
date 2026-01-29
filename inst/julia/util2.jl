## can maybe dlete more functions from this, chec this out later
#= CAN DELETE BELOW UP TO NEXT COMMENT IF THE CODE WORKS
const LIBNN = joinpath(@__DIR__, "libnn")  # i added this for r package 1/21/25 ??? do i need to change this in the util2-checkpoint file?????

# build nearest neighbor for observed location set S #
function BuildNN(coords, m, nThreads = 1.0, brute = 0.0)
    # coords: n by 2 array
    # m:      number of neighbors
    n = size(coords)[1]
    # coords_fit = [coords[1, :] coords[2, :]]
    nIndx = trunc(Cint,(1 + m) / 2 * m + (n - m - 1) * m);
    nnIndx = zeros(Cint, nIndx);
    nnDist = zeros(Float64, nIndx);
    ##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP)
    nnIndxLU = zeros(Cint, 2 * n);
    ##Note I need to convert n and m to pointers only because the R example of .C wants only pointers so that's how I worte the mkNNIndx and mkNNIndxTree0
    nPtr = Cint[n];
    mPtr = Cint[m];
    nThreadsPtr = Cint[nThreads];   

## i commented this out because i am changing it to the code below (needed it for r package) (replaced file path with libnn)
    #=if brute == 1.0
        ccall((:mkNNIndx, "../../utils/julia-R-nn-ccall2/libnn.dylib"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    else
        ccall((:mkNNIndxCB, "../../utils/julia-R-nn-ccall2/libnn.dylib"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    end=#

    if brute == 1.0
        ccall((:mkNNIndx, LIBNN), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    else
        ccall((:mkNNIndxCB, LIBNN), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    end
    
    nnIndx = nnIndx .+ 1 
    nnIndxLU = nnIndxLU[1:(n + 1)]  .+ 1
    nnIndxLU[n + 1] = nIndx + 1
    return (nnIndx = nnIndx, nnDist = nnDist, nnIndxLU = nnIndxLU[1:(n + 1)])
end  

=#




function getAD(coords, nnIndx, nnDist, nnIndxLU, ϕ, ν, A, D)
    # coords: n by 2 array,
    # nnIndx, nnDist, nnIndxLU output of the mkNNIndx function
    
    # return: A D
    # A = []
    # D = ones(size(coords)[1])
    
    # initialize memory
    NNdistM = []
    NNdistv = []
    n = length(nnIndxLU) - 1;
    if nnIndxLU[1] == nnIndxLU[2]
        D[1] = 1.0
        l = 2
    else 
        l = 1
    end
    for i in l:n
        if ν == 0.5
            NNdistM = cholesky(exp.(-ϕ .* pairwise(Euclidean(), 
                        coords[nnIndx[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)], :], 
                        dims = 1)))
            NNdistv = NNdistM.L \ (exp.(-ϕ .* nnDist[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)]))
            D[i] = 1.0 - NNdistv⋅NNdistv
            A[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)] = NNdistM.U \ NNdistv
        end
    end
    # return (A = A, D = D)
end

function getAD_collapse(coords, nnIndx, nnDist, nnIndxLU, ϕ, ν, α, A, D)
    # coords: n by 2 array,
    # nnIndx, nnDist, nnIndxLU output of the mkNNIndx function
    # ϕ, ν, α: covariance parameter set
    # hold: whether AD is for holded locations or not
    
    # return: A D
    # A = []
    # D = ones(size(coords)[1])
    
    # initialize memory
    NNdistM = []
    NNdistv = []
    n = length(nnIndxLU) - 1;
    if nnIndxLU[1] == nnIndxLU[2]
        D[1] = 1.0 / α
        l = 2
    else 
        l = 1
    end
    for i in l:n
        if ν == 0.5
            NNdistM = cholesky(exp.(-ϕ .* pairwise(Euclidean(), 
                        coords[nnIndx[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)], :], 
                        dims = 1)) + (1.0 / α - 1.0) * I)
            NNdistv = NNdistM.L \ (exp.(-ϕ .* nnDist[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)]))
            D[i] = 1.0 / α - NNdistv⋅NNdistv
            A[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)] = NNdistM.U \ NNdistv
        end
    end
    
    # return (A = A, D = D)
end


function kfoldperm(N,k)

    # k folder split of 1:N
    
    n,r = divrem(N,k)
    b = collect(1:n:N+1)
    for i in 1:length(b)
        b[i] += i > r ? r : i-1  
    end
    p = randperm(N)
    return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k]]
end

colnorm(A) = [norm(A[:,i], 2) for i=1:size(A,2)]

# Initialization #
using LinearAlgebra
using GeoStats
using Distances
function reorder_svd_by_spatial_range(Y_Xm::Matrix{Float64}, coords::Matrix{Float64}, K::Int, 
        n_s::Int = 1000, seed::Int = 123)
    # Set random seed
    Random.seed!(seed)

    n = size(Y_Xm, 1);
    n_sam = min(n_s, n);
    # Randomly sample n_sam rows
    sample_indices = randperm(n)[1:n_sam]
    Y_Xm_sample = Y_Xm[sample_indices, :]
    coords_sample = coords[sample_indices, :]

    # Compute randomized SVD
    U, S, Vt = rsvd(Y_Xm_sample, K);
    
    # Automatically compute maxlag based on coordinates
    max_dist = maximum(pairwise(Euclidean(), coords_sample, dims=1))
    
    # Compute ranges for each column of U
    ranges = Float64[]
    
    for k in 1:K
        # Create spatial data structure
        dat = georef((pc=U[:,k], x = coords_sample[:, 1], y = coords_sample[:, 2]), [:x, :y]);
        
        # Compute empirical variogram
        EV = EmpiricalVariogram(dat, :pc, 
            nlags=30,           # number of lags
            maxlag=0.5 * max_dist,         # maximum lag
            distance=Euclidean(), 
            estimator=:cressie  # default estimator
        )
        mod = GeoStatsFunctions.fit(ExponentialVariogram, EV);
        push!(ranges, range(mod).val);
    end
    
    # Sort indices based on ranges in descending order
    sorted_indices = sortperm(ranges, rev=true)
    
    # Reorder S and Vt
    S_reordered = S[sorted_indices]
    Vt_reordered = Vt[:, sorted_indices]
    
    println("Spatial ranges in descending order:")
    for (i, range) in enumerate(ranges[sorted_indices])
        println("PC$i: $range")
    end
    
    # Return Diagonal(S) * Vt'
    return Diagonal(S_reordered) * Vt_reordered' ./ sqrt(n_sam)
end


# diagonostic #
using DataFrames
using PrettyTables
using LinearAlgebra
using Statistics
function generate_comparison_table(
    F_mean_euclidean::Matrix{Float64}, 
    F_mean_projected::Matrix{Float64}, 
    F_proj::Matrix{Float64}, 
    digits::Int = 4
)
    # Validate input
    @assert size(F_mean_projected) == size(F_proj) == size(F_mean_euclidean) "Matrices must have same dimensions"
    
    n, K = size(F_mean_projected)
    
    # Compute metrics
    euclidean_distance = [
        round(norm(F_mean_projected[:, k] - F_proj[:, k]), digits=digits) 
        for k in 1:K
    ]
    
    # Sphere Variance: n-1 - ||mean||²
    sphere_variances = [
        round(n-1 - norm(F_mean_euclidean[:, k])^2, digits=digits) 
        for k in 1:K
    ]
    
    # Create DataFrame
    df = DataFrame(
        Metric = ["Euclidean Distance", "Sphere Var"]
    )
    
    # Add factor columns dynamically
    for k in 1:K
        df[!, "Factor $k"] = [euclidean_distance[k], sphere_variances[k]]
    end
    
    # Print table
    pretty_table(df)
    
    return df
end


