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
    # Check tensor shape
    two_b_geometry = size(tensor)
    D = length(two_b_geometry) ÷ 2

    @assert two_b_geometry[1:D] == V.geometry  "The tensor does not match the geometry of the lattice, geometry is $(V.geometry) vs tensor geometry $(two_b_geometry[1:D])"
    @assert two_b_geometry[1:D] == two_b_geometry[D+1:end] "The tensor does not satisfy two-body properties: second half of indices must match the first half"

    map_v_s = lattice.sites_v
    sites = collect(keys(map_v_s))

    Op = ZeroFockOperator()
    iszero(tensor) && return Op

    # Loop over all combinations of bra and ket sites
    for site_bra in sites
        for site_ket in sites
            ind = vcat(
                collect(map_v_s[site_bra]), 
                collect(map_v_s[site_ket])
                )
            coeff = tensor[ind...]
            if coeff != 0
                Op += FockOperator(
                    ((site_bra, true), (site_ket, false)),
                    coeff,
                    V
                )
            end
        end
    end

    return Op 
end


function four_body_Op(V::U1FockSpace, lattice::Lattice, tensor::AbstractArray{ComplexF64})
    # Check tensor shape
    four_b_geometry = size(tensor)
    D = length(four_b_geometry) ÷ 4

    @assert four_b_geometry[1:D] == V.geometry  "The tensor does not match the lattice geometry: $(V.geometry) vs $(four_b_geometry[1:D])"
    @assert four_b_geometry[1:D] == four_b_geometry[D+1:2*D] == four_b_geometry[2*D+1:3*D] == four_b_geometry[3*D+1:4*D]  "Tensor does not have consistent 4-body dimensions"

    map_v_s = lattice.sites_v
    sites = collect(keys(map_v_s))

    
    Op = ZeroFockOperator()
    iszero(tensor) && return Op
    # Loop over all combinations of 4 sites (bra1, bra2, ket1, ket2)
    for site1 in sites
        for site2 in sites
            for site3 in sites
                for site4 in sites
                    ind = vcat(
                        collect(map_v_s[site1]),  # bra1
                        collect(map_v_s[site2]),  # bra2
                        collect(map_v_s[site3]),  # ket1
                        collect(map_v_s[site4])   # ket2
                    )
                    coeff = tensor[ind...]
                    if coeff != 0
                        Op += FockOperator(
                            ((site1,true), (site2,true), (site3,false), (site4,false)),
                            coeff,
                            V
                        )
                    end
                end
            end
        end
    end

    return Op
end


function get_tensor_2body(Op::MultipleFockOperator, lattice::Lattice)
    #@warn "note that this function only takes into account terms of the form a†ᵢaⱼ and ignores all others"
    map_v_s = lattice.sites_v
    V = Op.terms[1].space
    geometry = V.geometry
    two_body_geometry = vcat(collect(geometry), collect(geometry)) |> Tuple
    tensor = zeros(ComplexF64, two_body_geometry...)

    for O in Op.terms
        if (length(O.product) == 2 ) & (O.product[1][2] & !O.product[2][2])
            s = map_v_s[O.product[1][1]]
            n = map_v_s[O.product[2][1]]
            ind = vcat(collect(s), collect(n))
            tensor[ind...] = O.coefficient 
        end
    end

    return tensor

end

function get_tensor_4body(Op::MultipleFockOperator, lattice::Lattice)
    #@warn "this function builds a tensor for 4-body terms a†_i a†_j a_k a_l"
    
    map_v_s = lattice.sites_v
    V = Op.terms[1].space
    geometry = V.geometry
    
    # Define 4-body tensor shape: (bra1, bra2, ket1, ket2)
    four_body_geometry = vcat(collect(geometry), collect(geometry), collect(geometry), collect(geometry)) |> Tuple
    tensor = zeros(ComplexF64, four_body_geometry...)

    for O in Op.terms
        # Select 4-body terms
        if length(O.product) == 4
            # Check operator structure: 2 daggers followed by 2 annihilation operators
            daggers = sum(p[2] for p in O.product[1:2])
            annih = sum(!p[2] for p in O.product[3:4])
            if (daggers == 2) & (annih == 2)
                # Map lattice sites to indices
                bra1 = map_v_s[O.product[1][1]]
                bra2 = map_v_s[O.product[2][1]]
                ket1 = map_v_s[O.product[3][1]]
                ket2 = map_v_s[O.product[4][1]]
                
                # Build tensor index
                ind = vcat(collect(bra1), collect(bra2), collect(ket1), collect(ket2))
                
                # Assign coefficient
                tensor[ind...] = O.coefficient
            end
        end
    end

    return tensor
end




function momentum_space_Op(Op::MultipleFockOperator, lattice::Lattice, dimensions::Tuple)
    V = Op.terms[1].space
    geometry = V.geometry
    D = length(geometry)

    dimensions_bra = dimensions
    dimensions_ket  = Tuple(collect(dimensions) .+ D)

    # --- 2-body ---
    real_tensor_2body = get_tensor_2body(Op, lattice)
    if iszero(real_tensor_2body)
        tensor_2body_m = zeros(ComplexF64,  nbody_geometry(geometry, 2))
    else
        tensor_2body_m = fft(real_tensor_2body, dimensions_bra)
        tensor_2body_m = ifft(tensor_2body_m, dimensions_ket) 
        tensor_2body_m = tensor_2body_m
    end

    # --- 4-body ---
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
        tensor_4body_m = tensor_4body_m
    end

    return two_body_Op(V, lattice, tensor_2body_m) + four_body_Op(V, lattice, tensor_4body_m)
end



end;