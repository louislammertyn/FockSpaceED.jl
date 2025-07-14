using FockSpace, JLD2

@time for N in 11:20, L in 1:5
    geometry = (L,)
    cutoff = N
    V = U1FockSpace(geometry, cutoff, N)
    basis(V)
end
