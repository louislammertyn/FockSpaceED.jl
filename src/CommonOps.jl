#####
#This file contains implementation of a selection of common Fock Operators.
#The first function calculates the general many-body wavefunction coefficient tensor as in the expression |ψ⟩ = ∑ Cₙₘₗ |n,m,l⟩
#Further impementations include:
#    - The 1-body density operator.
#    - The Entanglement Entropy.

begin
    
function MB_tensor(MBstate::MultipleFockState )
    s = MBstate.states[1]
    V = s.space
    modes = V.num_modes
    dims = ntuple(i-> (V.cutoff+1) , modes)
    C = zeros(ComplexF64, dims)
    for state in MBstate.states
        index = collect(state.occupations) .+1
        C[index...] = state.coefficient
    
    end
    return C
end

function Entanglement_Entropy(C::Array{ComplexF64, N}, cut::Int64) where N
    dims = size(C)
    d = dims[1]
    C_matrix = zeros(ComplexF64, d^cut, d^(N - cut))

    for i in CartesianIndices(C)
        row = 0
        for t in 0:(cut-1)
            ind = cut -t 
            row += (i[ind]-1) *  d^(t)
        end
        
        column = 0
        for t in 0:(N-cut -1)
            ind = N -t 
            column += (i[ind]-1) * d^(t)
        end
        
        C_matrix[row+1, column+1] =  C[i]
    end
    
    l=0
    # SVD and compute entropy
    _, S, _ = svd(C_matrix)
    p = S.^2 ./ sum(S.^2)  # Schmidt probabilities
    
    # Von Neumann entropy
    S_ent = -sum(p[p .> 0] .* log.(p[p .> 0])) # eps to avoid log(0)
    return S_ent, S
end

function density_onsite(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    matrix = zeros(ComplexF64, geometry)
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false)), 1. + 0im)
        matrix[s...] = state * (n * state)
    end
    return matrix
end

function density_flucs(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    matrix = zeros(ComplexF64, geometry)
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false), (sites[s], true), (sites[s], false)), 1. + 0im)
        matrix[s...] = state * (n * state) 
    end
    
    return matrix - (density_onsite(state, sites, geometry).^2)
end

function one_body_ρ(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    size_m = Tuple(vcat(collect(geometry), collect(geometry)))
    ρ = zeros(ComplexF64, size_m)
    for s1 in keys(sites), s2 in keys(sites)
        ind = vcat(collect(s1), collect(s2))
        Op = FockOperator(((sites[s1], true), (sites[s2], false)), 1. + 0im) 
        ρ[ind...] = state * (Op * state)
    end
    return ρ
end

function momentum_density(rho)
    
end

end