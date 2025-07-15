using Revise
using FockSpace, VectorInterface, KrylovKit

N = 25
geometry = (5,)
V = U1FockSpace(geometry, N,N)

basisFS(V)
s = fock_state(V, [1, 1, 3], 2. + 3im)
w = fock_state(V, [3, 0, 2], 1. + 4im)
scalartype(s)

scale(s, 2)
scale(ZeroFockState(),  2)

H = FockOperator(((1, true), (2, false) ), 1. +0im, V)
H +=  FockOperator(((2, true), (1, false) ), 1. +0im, V)
H+=  FockOperator(((2, true), (3, false) ), 1. +0im, V)
H+=  FockOperator(((3, true), (2, false) ), 1. +0im, V)
eigsolve(H, s, T=ComplexF64)
Tx = typeof(s)
Tfx = Core.Compiler.return_type(apply, Tuple{typeof(H),Tx})
T = Core.Compiler.return_type(dot, Tuple{Tx,Tfx})
@which dot