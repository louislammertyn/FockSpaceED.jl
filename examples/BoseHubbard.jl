using Revise
using FockSpace
using Test
using LinearAlgebra
using Plots
using ProgressMeter

J = 0.01
U = 1.
eigen([1 0; 3 2])

N = 8
geometry = (4,)
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

M = calculate_matrix_elements(states,H)

es, vs = eigen(M)
es
vs[:,1]
plot(1:length(vs[:,1]), abs2.(vs[:,1]))