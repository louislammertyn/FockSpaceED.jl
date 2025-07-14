using FockSpace

N = 4
geometry = (2,2, 2, 2)
latt = Lattice(geometry)

V = U1FockSpace(geometry, N,N)
s = fock_state(V, [1 1; 0 2], 1)

ms = MutableFockState(s)
a_j!(ms, 2)


