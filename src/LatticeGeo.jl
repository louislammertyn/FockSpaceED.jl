begin

abstract type AbstractLattice end

struct Lattice <:AbstractLattice
    D::Int
    sites::Dict{Tuple, Int}
    sites_v::Dict{Int,Tuple}
    NN::Dict{NTuple{Dim, Int}, Vector{NTuple{Dim,Int}}} where Dim

    function Lattice(geometry::NTuple{D,Int}; periodic::NTuple{D, Bool}=ntuple(i->false,D)) where D

        map_s_v = lattice_vectorisation_map(geometry)
        map_v_s = vector_to_lattice(map_s_v)
        NN = Lattice_NN(geometry; periodic)

        return new(D, map_s_v, map_v_s, NN)
    end 

end

############### Here we define the lattice functionalities for vectorisation of the lattice ###############

function vectorise_lattice(x::NTuple{D,Int}, L::NTuple{D,Int}) where D
    idx = 0
    stride = 1
    for i in 1:D
        idx += x[i] * stride
        stride *= L[i]
    end
    return idx +1
end

function lattice_vectorisation_map(geometry::NTuple{D,Int}) where D
    cartesian_to_vec = Dict{NTuple{D, Int}, Int}()
    for idx in Iterators.product(ntuple(i -> 0:geometry[i]-1, D)...)
        idx_s = ntuple(i-> idx[i]+1, D)
        cartesian_to_vec[idx_s] = vectorise_lattice(idx, geometry)
    end
    return cartesian_to_vec
end

vector_to_lattice(d::Dict) = Dict(v => k for (k, v) in d)

function Lattice_NN(geometry::NTuple{D, Int}; periodic::NTuple{D,Bool}=ntuple(i->false,D)) where D
    neighbours = Dict{NTuple{D, Int}, Vector{NTuple{D,Int}}}()
    ranges = (1:g for g in geometry)

    for site in Iterators.product(ranges...)
        site_tuple = Tuple(site)
        neighbours[site_tuple] = NTuple{D, Int}[]

        for d in 1:D
            for δ in (-1, 1)
                neighbor = collect(site)
                neighbor[d] += δ

                if 1 <= neighbor[d] <= geometry[d]
                    # within bounds
                    push!(neighbours[site_tuple], Tuple(neighbor))
                elseif periodic[d]
                    # periodic BCs
                    neighbor[d] = mod1(neighbor[d], geometry[d])
                    push!(neighbours[site_tuple], Tuple(neighbor))
                end
            end
        end
    end
    return neighbours
end


end;

lattice_vectorisation_map((4,4))