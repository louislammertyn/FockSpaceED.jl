begin

abstract type AbstractLattice end

struct Lattice <:AbstractLattice
    D::Int
    sites::Dict{Tuple, Int}
    sites_v::Dict{Int,Tuple}
    NN::Dict{NTuple{Dim, Int}, Vector{NTuple{Dim,Int}}} where Dim
    NN_v::Dict{Int,Vector{Int}}

    function Lattice(geometry::NTuple{D,Int}; periodic::NTuple{D, Bool}=ntuple(i->false,D), helical::NTuple{D, Bool}=ntuple(i->false,D)) where D
        map_s_v = lattice_vectorisation_map(geometry)
        map_v_s = vector_to_lattice(map_s_v)
        NN = (sum(helical) == false ? Lattice_NN(geometry; periodic) : Lattice_NN_h(geometry; periodic, helical))
        NN_v = vectorise_NN(NN, map_s_v)
        return new(D, map_s_v, map_v_s, NN, NN_v)
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

function Lattice_NN_h(geometry::NTuple{D, Int}; periodic::NTuple{D,Bool}=ntuple(i->false,D), helical::NTuple{D,Bool}=ntuple(i->false,D)) where D
    @assert D==2 "Code is implemented only for 2D setups with helical boundary conditions"
    for (ind,b) in enumerate(helical)
        if b
            @assert (!periodic[ind])  "The given boundary conditions are inconsistent (in at least one dimension periodic is false and helical true)"
        end
    end
    NN = Lattice_NN(geometry; periodic)

    for (ind, p) in enumerate(helical)
        !p && continue
        for site in keys(NN)
            ind_ = mod(ind,2)+1
            if (site[ind] == geometry[ind]) 
                site[ind_] == geometry[ind_] && continue
                println(site[ind]== geometry[ind])
                n = collect(site)
                n[ind] = 1
                n[ind_] += 1

                push!(NN[site], tuple(n...))
                push!(NN[tuple(n...)], site)
            end
        end
    end
    # Extract all sites
    sites = collect(keys(NN))
    x_sites = [s[1] for s in sites]
    y_sites = [s[2] for s in sites]

    # Start plotting
    pl = plot(; legend=false, xlim=(0, maximum(x_sites)+1), ylim=(0, maximum(y_sites)+1), aspect_ratio=1)
    
    # Draw lines for neighbors
    for (site, neighbors) in NN
        for n in neighbors
            # Plot a line between site and neighbor
            plot!(pl, [site[1], n[1]], [site[2], n[2]], color=:red, lw=1)
        end
    end
    scatter!(pl, x_sites, y_sites, color=:blue, markersize=6)
    display(pl)
    return NN        
end

function vectorise_NN(NN, map_s_v)
    NN_v = Dict()
    for s in keys(NN)
        sv = map_s_v[s]
        NN_v[sv] = []
        for n in NN[s]
            push!(NN_v[sv], map_s_v[n])
        end
    end
    return NN_v
end

end;

lattice_vectorisation_map((4,4))