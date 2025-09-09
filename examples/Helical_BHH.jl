using Revise
using FockSpace
using LinearAlgebra
using Plots
using ProgressMeter
using KrylovKit
using SparseArrays




N = 4
Lx = 4
Ly = 4

geometry = (Lx,Ly)
lattice = Lattice(geometry; periodic=(false,false), helical=(false,true));

lattice.NN
lattice.NN_v
lattice.sites
V = U1FockSpace(geometry, N, N)
states = all_states_U1_O(V)
#states = all_states_U1(V)


pl = plot(legend=false,ylims=[-.1,.1]);

for ϕ in LinRange(0,2*π, 5)

tx = exp(1im *ϕ)
ty=1 + 0im
U = 0

H = ZeroFockOperator()

for s in lattice.sites
    for n in lattice.NN[s[1]]
        if (n[1] > s[1][1])  & (n[2]==s[1][2])
            H += FockOperator(((s[2], true) , (lattice.sites[n], false)), tx, V)
        elseif ((n[2] > s[1][2]) & (n[1] == s[1][1])) || ((s[1][2]==geometry[2]) & (n[2] == 1) & (s[1][1]==n[1]-1))
            H += FockOperator(((s[2], true) , (lattice.sites[n], false)), ty, V)
        end
    end
end
H += dagger_FO(H)


M = calculate_matrix_elements_parallel(states,H)

es, vs = eigen(Hermitian(M))

i = findfirst(x->isapprox(x,0.; atol=01e-11), es)
gs = MultipleFockState(vs[:,i] .* states)
d = density_onsite(gs, lattice.sites, geometry)
display(heatmap(real.(d)))
scatter!(pl, ϕ.*ones(length(es)), es , marker=:o, color=:red, markersize=1, markerstrokecolor=:red)
end;

display(pl)

