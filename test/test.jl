using Revise
using FockSpace

N = 20
geometry = (5,)
latt = Lattice(geometry)

V = U1FockSpace(geometry, N, N)
Kin, Int = Bose_Hubbard_H(V, latt)
H = Kin + Int




states = basisFS(V)

@time M = calculate_matrix_elements_parallel(states, H);
@time M_ = calculate_matrix_elements_parallel_(states, H);
