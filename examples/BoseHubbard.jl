using Revise
using FockSpace
using Test
using LinearAlgebra
using Plots
using ProgressMeter

pd_gap = zeros(14,20)
pd = zeros(14,20)
@showprogress for n in 2:15, ju in 1:20
U = 1
J = (ju*0.02) * U


N = n
geometry = (3,)
D=length(geometry)

V = U1FockSpace(prod(geometry),N,N)
states = all_states_U1(geometry, V)

ind_v = lattice_vectorisation_map(geometry)
NN = Lattice_NN(geometry; periodic=(false,))

hoppings = zeros(ComplexF64, ntuple(i->geometry[mod(i,D)+1] , 2*D))
for site in keys(NN), n in NN[site]
    index = (site..., n...)
    hoppings[index...] = J
end 

H = ZeroFockOperator()

for site in keys(NN)
    for n in NN[site]
        index = (site..., n...)
        H += FockOperator(((ind_v[site], true), (ind_v[n], false)), hoppings[index...])
    end
end
for site in keys(NN)
    i = ind_v[site]
    H += FockOperator(((i, true ), (i,true), (i, false), (i, false)), U)
end

M = calculate_matrix_elements_parallel(states,H)

es, vs = eigen(Hermitian(M))
gs_coeff = vs[:,1]
pd_gap[(n-1),ju] =  (es[2]- es[1]) 
gs = create_MFS(states, gs_coeff)

ρ = zeros(3,3)

for i in 1:3, j in 1:3
    ρ_ij = FockOperator(((i, true), (j,false)), 1. +0im)
    ρ[i,j] = gs * (ρ_ij*gs)
end
es_rho, vs_rho = eigen(ρ)
pd[(n-1),ju] = es_rho[end] / n
end
heatmap(collect(1:20 ) .* 0.02, collect(2:15), pd, xlabel="J/U", ylabel="N", color=:viridis, title="Largest eigenvalue ρ")
heatmap(collect(1:20 ) .* 0.02, collect(2:15), pd_gap , xlabel="J/U", ylabel="N", color=:magma, title="Many body gap Δ")

