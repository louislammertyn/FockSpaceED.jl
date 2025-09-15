############################################################
# Fock Operator Utilities for Many-Body Quantum Systems
# 
# This file contains implementations of common Fock operators
# and utility functions to compute properties of many-body
# quantum states, such as:
#   - Full coefficient tensor representation of the state
#   - One-body and two-body density matrices
#   - On-site densities and fluctuations
#   - Entanglement entropy
#   - Mapping operators to momentum space
#   - Construction of common Hamiltonians (Bose-Hubbard)
############################################################

begin

############################################################
# Convert a MultipleFockState to a full many-body coefficient tensor
############################################################
"""
    MB_tensor(MBstate::MultipleFockState) -> Array{ComplexF64,N}

Given a `MultipleFockState` representing a many-body quantum state, 
returns the coefficient tensor `C` such that:

    |ψ⟩ = ∑ C[n,m,l,...] |n,m,l,...⟩

The dimensions of `C` are determined by the number of modes and the
cutoff in each mode.
"""
function MB_tensor(MBstate::MultipleFockState)
    s = MBstate.states[1]
    V = s.space
    modes = prod(V.geometry)
    dims = ntuple(i -> (V.cutoff + 1), modes)
    C = zeros(ComplexF64, dims)
    
    for state in MBstate.states
        index = collect(state.occupations) .+ 1
        C[index...] = state.coefficient
    end
    
    return C
end

############################################################
# Entanglement entropy via Schmidt decomposition
############################################################
"""
    Entanglement_Entropy(C::Array{ComplexF64,N}, cut::Int64) -> (S_ent, S)

Computes the von Neumann entanglement entropy for a bipartition
of the system after reshaping the coefficient tensor `C`.

- `cut`: number of modes in subsystem A.
- Returns:
    - `S_ent`: entanglement entropy
    - `S`: singular values (Schmidt coefficients)
"""
function Entanglement_Entropy(C::Array{ComplexF64,N}, cut::Int64) where N
    dims = size(C)
    d = dims[1]
    C_matrix = zeros(ComplexF64, d^cut, d^(N - cut))

    # Map N-dimensional tensor indices to 2D matrix
    for i in CartesianIndices(C)
        row = 0
        for t in 0:(cut-1)
            ind = cut - t
            row += (i[ind]-1) * d^t
        end

        column = 0
        for t in 0:(N - cut - 1)
            ind = N - t
            column += (i[ind]-1) * d^t
        end

        C_matrix[row + 1, column + 1] = C[i]
    end

    # Compute singular values and probabilities
    _, S, _ = svd(C_matrix)
    p = S.^2 ./ sum(S.^2)  # Schmidt probabilities

    # Von Neumann entropy
    S_ent = -sum(p[p .> 0] .* log.(p[p .> 0]))
    return S_ent, S
end

############################################################
# On-site density operator
############################################################
"""
    density_onsite(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int}) -> Array{ComplexF64,D}

Computes the expectation value ⟨n_i⟩ on each lattice site for a given Fock state.
"""
function density_onsite(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    matrix = zeros(ComplexF64, geometry)
    V = typeof(state) == MultipleFockState ? state.states[1].space : state.space
    
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false)), 1. + 0im, V)
        matrix[s...] = state * (n * state)
    end
    return matrix
end

function center_of_mass(densities::AbstractArray)
    geometry = size(densities)
    CoM = zeros(length(geometry))
    for (i,L) in enumerate(geometry)
        shape = ntuple(d -> d==i ? L : 1, ndims(densities))
        w = reshape(collect(1:L), shape)
        CoM[i] = sum(densities .* w) / sum(densities)
    end
    return CoM
end

############################################################
# On-site density fluctuations
############################################################
"""
    density_flucs(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int}) -> Array{ComplexF64,D}

Computes the variance ⟨n_i^2⟩ - ⟨n_i⟩^2 for each site.
"""
function density_flucs(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    matrix = zeros(ComplexF64, geometry)
    for s in keys(sites)
        n = FockOperator(((sites[s], true), (sites[s], false),
                          (sites[s], true), (sites[s], false)), 1. + 0im, state.space)
        matrix[s...] = state * (n * state)
    end
    return matrix - density_onsite(state, sites, geometry).^2
end

############################################################
# One-body density matrix
############################################################
"""
    one_body_ρ(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int}) -> Array{ComplexF64,2D}

Computes the one-body density matrix ρ_{ij} = ⟨a_i^† a_j⟩.
"""
function one_body_ρ(state::AbstractFockState, sites::Dict, geometry::NTuple{D, Int64}) where D
    V = typeof(state) == MultipleFockState ? state.states[1].space : state.space
    size_m = Tuple(vcat(collect(geometry), collect(geometry)))
    ρ = zeros(ComplexF64, size_m)

    for s1 in keys(sites), s2 in keys(sites)
        ind = vcat(collect(s1), collect(s2))
        Op = FockOperator(((sites[s1], true), (sites[s2], false)), 1. + 0im, V)
        ρ[ind...] = state * (Op * state)
    end

    return ρ
end

############################################################
# Hamiltonians: Bose-Hubbard
############################################################
"""
    Bose_Hubbard_H(V::U1FockSpace, lattice::Lattice, J::Number=1., U::Number=1.) -> (Kin, Int)

Constructs the kinetic and interaction parts of the Bose-Hubbard Hamiltonian:

- `Kin`: hopping term H_J
- `Int`: on-site interaction term H_U
"""
function Bose_Hubbard_H(V::U1FockSpace, lattice::Lattice, J::Number=1., U::Number=1.)
    N = V.particle_number
    D = length(V.geometry)
    sites = lattice.sites
    NN = lattice.NN

    # Precompute hopping tensor
    hoppings = zeros(ComplexF64, ntuple(i -> V.geometry[mod(i,D)+1], 2*D))
    for site in keys(NN), n in NN[site]
        index = (site..., n...)
        hoppings[index...] = J
    end 

    Kin = ZeroFockOperator()
    Int = ZeroFockOperator()

    # Kinetic term
    for site in keys(NN), n in NN[site]
        index = (site..., n...)
        Kin += FockOperator(((sites[site], true), (sites[n], false)), hoppings[index...], V)
    end

    # Interaction term
    for site in keys(NN)
        i = sites[site]
        Int += FockOperator(((i, true ), (i,true), (i, false), (i, false)), U, V)
    end

    return Kin, Int
end

############################################################
# Two-body operator from tensor
############################################################
"""
    two_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64}) -> MultipleFockOperator

Constructs a two-body Fock operator from a 2-body tensor.
"""
function two_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64})
    two_b_geometry = size(tensor)
    D = length(two_b_geometry) ÷ 2

    @assert two_b_geometry[1:D] == V.geometry  "Tensor geometry mismatch"
    @assert two_b_geometry[1:D] == two_b_geometry[D+1:end] "Tensor is not 2-body symmetric"

    map_v_s = lattice.sites_v
    sites = collect(keys(map_v_s))
    Op = ZeroFockOperator()
    iszero(tensor) && return Op

    for site_bra in sites, site_ket in sites
        ind = vcat(collect(map_v_s[site_bra]), collect(map_v_s[site_ket]))
        coeff = tensor[ind...]
        if coeff != 0
            Op += FockOperator(((site_bra, true), (site_ket, false)), coeff, V)
        end
    end

    return Op
end

############################################################
# Four-body operator from tensor
############################################################
"""
    four_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64}) -> MultipleFockOperator

Constructs a four-body Fock operator from a 4-body tensor.
"""
function four_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64})
    four_b_geometry = size(tensor)
    D = length(four_b_geometry) ÷ 4

    @assert four_b_geometry[1:D] == V.geometry  "Tensor geometry mismatch"
    @assert four_b_geometry[1:D] == four_b_geometry[D+1:2*D] == four_b_geometry[2*D+1:3*D] == four_b_geometry[3*D+1:4*D]

    map_v_s = lattice.sites_v
    sites = collect(keys(map_v_s))
    Op = ZeroFockOperator()
    iszero(tensor) && return Op

    for site1 in sites, site2 in sites, site3 in sites, site4 in sites
        ind = vcat(collect(map_v_s[site1]), collect(map_v_s[site2]),
                   collect(map_v_s[site3]), collect(map_v_s[site4]))
        coeff = tensor[ind...]
        if coeff != 0
            Op += FockOperator(((site1,true),(site2,true),(site3,false),(site4,false)), coeff, V)
        end
    end

    return Op
end

############################################################
# Tensor extraction and momentum-space operators
############################################################

############################################################
# Extract 2-body tensor from a MultipleFockOperator
############################################################
"""
    get_tensor_2body(Op::MultipleFockOperator, lattice::Lattice) -> Array{ComplexF64,2D}

Constructs a 2-body tensor representation of terms of the form a†_i a_j
from a `MultipleFockOperator`. Ignores other types of terms.

Arguments:
- `Op`: MultipleFockOperator containing operator terms
- `lattice`: Lattice object mapping sites to indices

Returns:
- `tensor`: 2D (or 2*D-dimensional) array of complex coefficients
"""
function get_tensor_2body(Op::MultipleFockOperator, lattice::Lattice)
    map_v_s = lattice.sites_v
    V = Op.terms[1].space
    geometry = V.geometry
    two_body_geometry = Tuple(vcat(collect(geometry), collect(geometry)))
    tensor = zeros(ComplexF64, two_body_geometry...)

    # Loop over all terms
    for O in Op.terms
        # Select terms with exactly one creation and one annihilation operator
        if (length(O.product) == 2) & (O.product[1][2] & !O.product[2][2])
            s = map_v_s[O.product[1][1]]  # creation site
            n = map_v_s[O.product[2][1]]  # annihilation site
            ind = vcat(collect(s), collect(n))
            tensor[ind...] = O.coefficient
        end
    end

    return tensor
end

############################################################
# Extract 4-body tensor from a MultipleFockOperator
############################################################
"""
    get_tensor_4body(Op::MultipleFockOperator, lattice::Lattice) -> Array{ComplexF64,4D}

Constructs a 4-body tensor representation of terms of the form
a†_i a†_j a_k a_l from a `MultipleFockOperator`.

Arguments:
- `Op`: MultipleFockOperator containing operator terms
- `lattice`: Lattice object mapping sites to indices

Returns:
- `tensor`: 4D (or 4*D-dimensional) array of complex coefficients
"""
function get_tensor_4body(Op::MultipleFockOperator, lattice::Lattice)
    map_v_s = lattice.sites_v
    V = Op.terms[1].space
    geometry = V.geometry
    four_body_geometry = Tuple(vcat(collect(geometry), collect(geometry),
                                    collect(geometry), collect(geometry)))
    tensor = zeros(ComplexF64, four_body_geometry...)

    # Loop over all operator terms
    for O in Op.terms
        if length(O.product) == 4
            # Count creation and annihilation operators
            daggers = sum(p[2] for p in O.product[1:2])
            annih = sum(!p[2] for p in O.product[3:4])
            
            if (daggers == 2) & (annih == 2)
                # Map lattice sites to indices
                bra1 = map_v_s[O.product[1][1]]
                bra2 = map_v_s[O.product[2][1]]
                ket1 = map_v_s[O.product[3][1]]
                ket2 = map_v_s[O.product[4][1]]

                # Build tensor index and assign coefficient
                ind = vcat(collect(bra1), collect(bra2), collect(ket1), collect(ket2))
                tensor[ind...] = O.coefficient
            end
        end
    end

    return tensor
end

############################################################
# Map Fock operator to momentum space
############################################################
"""
    momentum_space_Op(Op::MultipleFockOperator, lattice::Lattice, dimensions::Tuple) -> MultipleFockOperator

Transforms a `MultipleFockOperator` to momentum space using FFTs.

Arguments:
- `Op`: MultipleFockOperator to transform
- `lattice`: Lattice object
- `dimensions`: tuple of FFT dimensions for each spatial direction

Returns:
- `Op_momentum`: MultipleFockOperator in momentum space
"""
function momentum_space_Op(Op::MultipleFockOperator, lattice::Lattice, dimensions::Tuple)
    V = Op.terms[1].space
    geometry = V.geometry
    D = length(geometry)
    dimensions_bra = dimensions
    dimensions_ket = Tuple(collect(dimensions) .+ D)

    # --- 2-body transformation ---
    real_tensor_2body = get_tensor_2body(Op, lattice)
    if iszero(real_tensor_2body)
        tensor_2body_m = zeros(ComplexF64, nbody_geometry(geometry, 2))
    else
        tensor_2body_m = fft(real_tensor_2body, dimensions_bra)
        tensor_2body_m = ifft(tensor_2body_m, dimensions_ket)
    end

    # --- 4-body transformation ---
    real_tensor_4body = get_tensor_4body(Op, lattice)
    if iszero(real_tensor_4body)
        tensor_4body_m = zeros(ComplexF64, nbody_geometry(geometry, 4))
    else
        bra_dims_4body  = collect(dimensions)
        bra_dims_4body2 = collect(dimensions) .+ D
        ket_dims_4body  = collect(dimensions) .+ 2*D
        ket_dims_4body2 = collect(dimensions) .+ 3*D

        tensor_4body_m = fft(real_tensor_4body, bra_dims_4body)
        tensor_4body_m = fft(tensor_4body_m, bra_dims_4body2)
        tensor_4body_m = ifft(tensor_4body_m, ket_dims_4body)
        tensor_4body_m = ifft(tensor_4body_m, ket_dims_4body2)
        tensor_4body_m .= tensor_4body_m
    end

    # Construct momentum-space operator
    return two_body_Op(V, lattice, tensor_2body_m) + four_body_Op(V, lattice, tensor_4body_m)
end


############################################################
# End of Fock operator utilities
############################################################

end;
