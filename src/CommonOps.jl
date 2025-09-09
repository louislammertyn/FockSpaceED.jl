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
    modes = prod(V.geometry)
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
    if typeof(state) == MultipleFockState
        V = state.states[1].space
    else
        V = state.space
    end
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false)), 1. + 0im, V)
        matrix[s...] = state * (n * state)
    end
    return matrix
end

function density_flucs(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    matrix = zeros(ComplexF64, geometry)
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false), (sites[s], true), (sites[s], false)), 1. + 0im, state.space)
        matrix[s...] = state * (n * state) 
    end
    
    return matrix - (density_onsite(state, sites, geometry).^2)
end

function one_body_ρ(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    typeof(state) == MultipleFockState ? V = state.states[1].space : V = state.space
    size_m = (Tuple(vcat(collect(geometry), collect(geometry))))
    ρ = zeros(ComplexF64, size_m)
    for s1 in keys(sites), s2 in keys(sites)
        ind = vcat(collect(s1), collect(s2))
        Op = FockOperator(((sites[s1], true), (sites[s2], false)), 1. + 0im, V) 
        ρ[ind...] = state * (Op * state)
    end
    return ρ
end

function momentum_density(rho)
    
end

############## Common Hamiltonians ##############

function Bose_Hubbard_H(V::U1FockSpace, lattice::Lattice, J::Number=1., U::Number=1.)
    N = V.particle_number
    D = length(V.geometry)

    sites = lattice.sites
    NN = lattice.NN

    hoppings = zeros(ComplexF64, ntuple(i->V.geometry[mod(i,D)+1] , 2*D))
    for site in keys(NN), n in NN[site]
        index = (site..., n...)
        hoppings[index...] = J
    end 

    Kin = ZeroFockOperator()
    Int = ZeroFockOperator()

    for site in keys(NN)
        for n in NN[site]
            index = (site..., n...)
            Kin += FockOperator(((sites[site], true), (sites[n], false)), hoppings[index...], V)
        end
    end
    for site in keys(NN)
        i = sites[site]
        Int += FockOperator(((i, true ), (i,true), (i, false), (i, false)), U, V)
    end
    return Kin, Int
end


function two_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64})
    two_b_geometry = size(tensor)
    D = length(two_b_geometry) / 2

    @assert two_b_geometry[1:D] == V.geometry  "The tensor does not match the geometry of the lattice, geometry is $(V.geometry) and tensor geometry is $(size(tensor)[1:D])" 
    @assert two_b_geometry[1:D] == two_b_geometry[D+1:end] "The tensor does not satisfy the properties of a two body tensor. Please check if the second half of tensori ndices has the same size as the first half."

    map_v_s = lattice.sites_v
    NN_v = lattice.NN_v
    Op = ZeroFockOperator()

    for site in keys(NN_v)
        for n in NN_v[site]
            Ind = vcat(collect(map_v_s[site]), collect(map_v_s[n])) 
            Op += FockOperator(((site, true), (n, false)), tensor[Ind...], V)
        end
    end

    return Op

end

function get_tensor_two_body(Op::MultipleFockOperator, lattice::Lattice)
    map_v_s = lattice.sites_v
    V = Op.terms[1].V
    geometry = V.geometry
    two_body_geometry = vcat(collect(geometry), collect(geometry)) |> Tuple
    tensor = zeros(ComplexF64, two_body_geometry)

    for O in Op.terms
        if length(O.product) == 2
            s = map_v_s[O.product[1][1]]
            n = map_v_s[O.product[2][1]]
            ind = vcat(collect(s), collect(n))
            tensor[ind...] = O.coefficient 
        end
    end
end


delta(i::Int, j::Int) = (i==j ? 1 : 0)
end;