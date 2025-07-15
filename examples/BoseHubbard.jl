using Revise
using FockSpace
using Test
using LinearAlgebra
using Plots
using ProgressMeter
using KrylovKit


############ Phase Diagram ####################


Nrange = 5:7
jurange = 1:2
pd_gap = zeros(Nrange[end],jurange[end])
pd = zeros(Nrange[end],jurange[end])
#@showprogress 
@time for n in Nrange, ju in jurange
U = 1
J = (ju*0.02) * U
N = n+1
E_scale = J+ U*N

geometry = (5,)
D=length(geometry)

V = U1FockSpace(geometry,N,N)
states = basisFS(V)

latt = Lattice(geometry)
ind_v = latt.sites
NN = latt.NN

hoppings = zeros(ComplexF64, ntuple(i->geometry[mod(i,D)+1] , 2*D))
for site in keys(NN), n in NN[site]
    index = (site..., n...)
    hoppings[index...] = J
end 

H = ZeroFockOperator()

for site in keys(NN)
    for n in NN[site]
        index = (site..., n...)
        H += FockOperator(((ind_v[site], true), (ind_v[n], false)), hoppings[index...], V)
    end
end
for site in keys(NN)
    i = ind_v[site]
    H += FockOperator(((i, true ), (i,true), (i, false), (i, false)), U, V)
end

M = calculate_matrix_elements_parallel(states,H)

@time es, vs = eigen(Hermitian(M))

gs_coeff = vs[:,1]
pd_gap[(n),ju] =  (es[2]- es[1]) / E_scale
gs = create_MFS( gs_coeff, states)

ρ = zeros(3,3)

for i in 1:3, j in 1:3
    ρ_ij = FockOperator(((i, true), (j,false)), 1. +0im, V)
    ρ[i,j] = gs * (ρ_ij*gs)
end
es_rho, vs_rho = eigen(ρ)
pd[(n),ju] = es_rho[end] / N
end;

heatmap(collect(jurange ) .* 0.02, collect(Nrange) .+ 1, pd, xlabel="J/U", ylabel="N", color=:viridis, title="Largest eigenvalue ρ")
heatmap(collect(jurange ) .* 0.02, collect(Nrange) .+ 1, pd_gap , xlabel="J/U", ylabel="N", color=:magma, title="Many body gap Δ")




##### 2D with Krylov methods


U = 1
J = .1* U


N = 5
geometry = (5,)


V = U1FockSpace(geometry,N,N)
states = basis(V)

latt = Lattice(geometry)
sites = latt.sites
NN = latt.NN

hoppings = zeros(ComplexF64, ntuple(i->geometry[mod(i,D)+1] , 2*D))
for site in keys(NN), n in NN[site]
    index = (site..., n...)
    hoppings[index...] = J
end 

H = ZeroFockOperator()

for site in keys(NN)
    for n in NN[site]
        index = (site..., n...)
        H += FockOperator(((sites[site], true), (sites[n], false)), hoppings[index...])
    end
end
for site in keys(NN)
    i = sites[site]
    H += FockOperator(((i, true ), (i,true), (i, false), (i, false)), U)
end

M = calculate_matrix_elements_parallel(states,H)

es, vs = eigen(Hermitian(M))
gs = create_MFS(vs[:,1], states)

density_onsite(gs, sites, geometry)
es, _ = eigen(one_body_ρ(gs, sites, geometry))
density_flucs(gs, sites, geometry)
rho = one_body_ρ(gs, sites, geometry)