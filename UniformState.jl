md"""
This file introduces some custom types that make the vectorspace structure of the Fock space compatible
with KrylovKit.

The main issue is that the AbstractFockStates have several subtypes and outcomes of operations of operators on
Fock states can result in several subtypes. To fix this problem we will make a single vectorlike object that encodes
the Fock space structure in a uniform type stable way.

"""
using LinearAlgebra, VectorInterface, FockSpace

struct U1FockVectorSpace <: AbstractFockSpace
    basis::NTuple{N, AbstractFockState} where N
end
function U1FockVectorSpace(states::Vector{AbstractFockState})
    return U1FockVectorSpace(ntuple(i -> states[i], length(states)))
end


mutable struct FockVector <: AbstractVector{ComplexF64}
    coefficients::Vector{ComplexF64}
    V::U1FockVectorSpace
end

FockVector(V::U1FockVectorSpace) = FockVector(zeros(ComplexF64, length(V.basis)), V)


# -- Required AbstractVector Interface --
Base.zero(v::FockVector) = v.coefficients .* 0
Base.size(v::FockVector) = (length(v.coefficients),)
Base.length(v::FockVector) = length(v.coefficients)
Base.IndexStyle(::Type{<:FockVector}) = IndexLinear()

Base.getindex(v::FockVector, i::Int) = v.coefficients[i]
Base.setindex!(v::FockVector, val, i::Int) = (v.coefficients[i] = val)

Base.eltype(::Type{FockVector}) = ComplexF64

# -- Optional: iteration, printing, dot product, etc. --

Base.iterate(v::FockVector, state=1) =
    state > length(v) ? nothing : (v[state], state + 1)

function Base.show(io::IO, v::FockVector)
    println(io, "FockVector of dimension $(length(v))")
    for (i, coeff) in enumerate(v)
        bs = v.V.basis[i]
        println(io, "  $coeff × |$(bs.occupations)⟩")
    end
end

# -- Optional utility constructor --

function FockVector(V::U1FockVectorSpace; coeffs=nothing)
    N = length(V.basis)
    if isnothing(coeffs)
        coeffs = zeros(ComplexF64, N)
    elseif length(coeffs) != N
        error("Length of coefficients must match dimension of basis.")
    end
    return FockVector(coeffs, V)
end

using VectorInterface
V_ = U1FockSpace((5,), 5,5)
V = U1FockVectorSpace(basisFS(V_))

a = FockVector(V)

function basis_v(i::Int64, V::U1FockSpace)
    l = prod(V.geometry)
    v = zeros(Int64, l)
end


function bounded_base_sums_parallel(base::Int, total::Int, maxlen::Int)
    max_i = base^maxlen  # this guarantees all digit sequences up to maxlen
    thread_results = [Vector{Vector{Int}}() for _ in 1:Threads.nthreads()]
    
    Threads.@threads for i in 0:max_i-1
        n = digits(i, base=base)
        if sum(n) == total && length(n) <= maxlen
            push!(thread_results[Threads.threadid()], n)
        end
    end

    raw = reduce(vcat, thread_results)
    padded = raw .|> x -> vcat(zeros(Int, maxlen - length(x)), x)
    sorted = sort(padded, by= x->evalpoly(6, x), rev=true)

    return sorted
end

N = 10
L = 10
@time b = bounded_base_sums_deterministic(N+1, N,L)

Threads.nthreads()
@time filter(i->(sum(i)==5), b)
maximum(length, b)

# Define VectorInterface.jl necessities
VectorInterface.scalartype(s::AbstractFockState) = ComplexF64

VectorInterface.zerovector(s::FockVector) = zero(s)
VectorInterface.zerovector!(s::MutableFockState) = 0. * s
VectorInterface.zerovector!!(s::AbstractFockState) = isa(s, MutableFockState) ? VectorInterface.zerovector!(s) : VectorInterface.zerovector(s)

VectorInterface.scale(s::AbstractFockState, α::Number) = α * s
VectorInterface.scale!(s::MutableFockState, α::Number) = mul_Mutable!(α, s)
VectorInterface.scale!!(s::AbstractFockState, α::Number) = isa(s, MutableFockState) ? VectorInterface.scale!(s, α) : VectorInterface.scale(s, α)

VectorInterface.add(s::AbstractFockState, w::AbstractFockState, α::Number=1. +0im,  β::Number=1. +0im) = α * s + β * w
VectorInterface.add!(s::AbstractFockState, w::AbstractFockState, α::Number=1. +0im,  β::Number=1. +0im) = @assert false "In-place addition is not defined for MutableFockState"
VectorInterface.add!!(s::AbstractFockState, w::AbstractFockState, α::Number=1. +0im,  β::Number=1. +0im) = isa(s, MutableFockState) ? VectorInterface.add!(s, α) : VectorInterface.add(s, α)

VectorInterface.inner(s1::AbstractFockState, s2::AbstractFockState) = s1 * s2
VectorInterface.norm(s::AbstractFockState) = sqrt(s * s)
LinearAlgebra.dot(a::FockVector, b::FockVector) = a.coefficients' * b.coefficients 
