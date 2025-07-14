
############### In this part we define the operators acting on a Fock space ########
# Each FockOperator consists of a product of creation and annihilation operators encoded as a tuple of tuples. 
# Each tuple corresponds to an (i,a/ad) with position i and {ad <-> true, a <-> false}
# The order of the tuple corresponds to the way we read the states on paper i.e. THE LAST TUPLE WILL BE THE ONE TO ACT FIRST ON THE STATE
# There is also a coefficient for the term
# Example:  θ = c1 * ad_1 * ad_2 * a_3 * a_2 * ad_6    <->     θ = FockOperator[( (1, true), (2, true), (3, false), (2, false), (6,true) ), c1]
# MultipleFockOperator represents sums of FockOperator types
begin

abstract type AbstractFockOperator end

struct FockOperator <: AbstractFockOperator
    product::NTuple{N, Tuple{Int,Bool}} where N 
    coefficient::ComplexF64
    space::AbstractFockSpace
end

mutable struct MultipleFockOperator <: AbstractFockOperator
    terms::Vector{FockOperator}
end

# Neutral element structure
struct ZeroFockOperator <: AbstractFockOperator
end

# empty product means identity

function identity_fockoperator(V::AbstractFockSpace, c::ComplexF64=1.0)
    FockOperator(NTuple{0, Tuple{Int,Bool}}() , c, V) 
end

Base.:+(z::ZeroFockOperator, s::FockOperator) = s
Base.:+(s::FockOperator, z::ZeroFockOperator) = s
Base.:+(z::ZeroFockOperator, ms::MultipleFockOperator) = ms
Base.:+(ms::MultipleFockOperator, z::ZeroFockOperator) = ms
Base.:-(z::ZeroFockOperator, s::FockOperator) = -1 *s
Base.:-(s::FockOperator, z::ZeroFockOperator) = s
Base.:-(z::ZeroFockOperator, ms::MultipleFockOperator) =-1* ms
Base.:-(ms::MultipleFockOperator, z::ZeroFockOperator) = ms
Base.:*(c::Number, z::ZeroFockOperator) = z
Base.:*(z::ZeroFockOperator, c::Number) = z
Base.:*(z1::ZeroFockOperator, z2::ZeroFockOperator) = 0
Base.:*(z::ZeroFockOperator, s::FockOperator) = 0
Base.:*(s::FockOperator, z::ZeroFockOperator) = 0

######## Pretty printing #########
function Base.show(io::IO, op::FockOperator)
    str = op.coefficient == 1 + 0im ? "" : string("($(op.coefficient))", " ⋅ ")

    # Build readable string
    for (site, is_creation) in op.product
        str *= is_creation ? "cr($site)" : "ann($site)"
        str *= " "
    end

    print(io, strip(str))
end

function Base.show(io::IO, mop::MultipleFockOperator)
    if isempty(mop.terms)
        print(io, "0")
        return
    end

    for (i, term) in enumerate(mop.terms)
        if i > 1
            print(io, " + ")
        end
        print(io, term)
    end
end

########## Basic operations ##########
Base.size(Op::FockOperator) = (prod(Op.space.geometry), prod(Op.space.geometry))

Base.size(Op::MultipleFockOperator) = size(Op.terms[1])

Base.eltype(Op::AbstractFockOperator) = ComplexF64

Base.:+(op1::FockOperator, op2::FockOperator) =
    op1.product == op2.product ? cleanup_FO(FockOperator(op1.product, op1.coefficient + op2.coefficient, op1.space)) : MultipleFockOperator([op1, op2]);

Base.:-(op1::FockOperator, op2::FockOperator) = op1 + FockOperator(op2.product, -op2.coefficient, op2.space)

function Base.:+(op::FockOperator, mop::MultipleFockOperator)
    new_terms = copy(mop.terms)
    matched = false
    for i in eachindex(new_terms)
        if new_terms[i].product == op.product
            new_terms[i] = FockOperator(op.product, new_terms[i].coefficient + op.coefficient, op.space)
            matched = true
            break
        end
    end
    if !matched
        push!(new_terms, op)
    end
    return cleanup_FO(MultipleFockOperator(new_terms))
end

Base.:+(mop::MultipleFockOperator, op::FockOperator) = op + mop
Base.:-(mop::MultipleFockOperator, op::FockOperator) = mop + FockOperator(op.product, -op.coefficient, op.space)
Base.:-(op::FockOperator, mop::MultipleFockOperator) = (-1) * mop + op

function Base.:+(mop1::MultipleFockOperator, mop2::MultipleFockOperator)
    result = copy(mop2)
    for t in mop1.terms
        result = result + t
    end
    return cleanup_FO(result)
end

Base.:-(mop1::MultipleFockOperator, mop2::MultipleFockOperator) = mop1 + (-1) * mop2

Base.:*(c::Number, op::FockOperator) = FockOperator(op.product, c * op.coefficient, op.space)
Base.:*(op::FockOperator, c::Number) = c * op

function Base.:*(c::Number, mop::MultipleFockOperator)
    new_terms = [c * t for t in mop.terms]
    return cleanup_FO(MultipleFockOperator(new_terms))
end

Base.:*(mop::MultipleFockOperator, c::Number) = c * mop

# Multiplying operators
function Base.:*(Op1::FockOperator, Op2::FockOperator)
    factors1 = collect(Op1.product)
    factors2 = collect(Op2.product)
    new_factor = vcat(factors1, factors2)
    return FockOperator(Tuple(new_factor), Op1.coefficient * Op2.coefficient, Op1.space)
end

function Base.:*(MOp::MultipleFockOperator, Op::FockOperator)
    terms = Vector{FockOperator}()
    for O in MOp.terms
        push!(terms, O * Op)
    end
    return MultipleFockOperator(terms)
end

function Base.:*(Op::FockOperator, MOp::MultipleFockOperator)
    terms = Vector{FockOperator}()
    for O in MOp.terms
        push!(terms, Op * O)
    end
    return MultipleFockOperator(terms)
end

function Base.:*(MOp1::MultipleFockOperator, MOp2::MultipleFockOperator)
    terms = Vector{FockOperator}()
    for O1 in MOp1.terms, O2 in MOp2.terms
        push!(terms, O1 * O2)
    end
    return MultipleFockOperator(terms)
end

########## Utilities ##########
function Base.copy(op::FockOperator)
    return FockOperator(op.product, op.coefficient, op.space)
end

Base.copy(mop::MultipleFockOperator) = MultipleFockOperator(copy(mop.terms))

function cleanup_FO(op::FockOperator)
    return op.coefficient==0. ? ZeroFockOperator() : op
end
function cleanup_FO(mop::MultipleFockOperator)
    new_terms = filter(t -> t.coefficient != 0, mop.terms)
    if length(new_terms) == 0
        return ZeroFockOperator()
    elseif length(new_terms) == 1
        return new_terms[1]
    else
        return MultipleFockOperator(new_terms)
    end
end

function cleanup_FO(mop::ZeroFockOperator)
    return mop
end

function dagger_FO(Op::FockOperator)
    c_dag = Op.coefficient'
    new_terms = []
    for o in reverse(Op.product)
        tup = (o[1],!o[2])
        push!(new_terms, tup)
    end
    return FockOperator(Tuple(new_terms), c_dag, Op.space)
end

function dagger_FO(Ops::MultipleFockOperator)
    new_ops = []
    for op in Ops.terms
        push!(new_ops, dagger_FO(op))
    end
    return MultipleFockOperator(new_ops)
end

function Base.:*(Op::FockOperator, ket::AbstractFockState)
    for factor in reverse(Op.product)
        if factor[2]
            ket = ad_j(ket, factor[1])
        else
            ket = a_j(ket, factor[1])
        end
    end 
    return Op.coefficient * ket
end

function Base.:*(Ops::MultipleFockOperator, ket::AbstractFockState)
    new_ket = ZeroFockState()

    for Op in Ops.terms
        new_ket = new_ket + (Op * ket)
        checkU1(new_ket)
    end

    return new_ket
end 

function apply!(Op::FockOperator, ket::MutableFockState)
    for (site, ladder) in reverse(Op.product)
        if ladder
            ad_j!(ket, site)
        else
            a_j!(ket, site)
        end
    end
    mul_Mutable!(Op.coefficient, ket)
end    


function calculate_matrix_elements(states::Vector{AbstractFockState}, Ops::MultipleFockOperator)
    Op_matrix = zeros(ComplexF64, length(states), length(states))
    tmp = MutableFockState(states[1])
    
    for (i,bra) in enumerate(states), (j,ket) in enumerate(states)  
        total_ij = 0.0 + 0im
        for Op in Ops.terms
            reset2!(tmp, ket.occupations, ket.coefficient)                 
            apply!(Op, tmp)   
            if tuple_vector_equal(bra.occupations, tmp.occupations)
                total_ij += bra.coefficient' * tmp.coefficient
            end
        end
        if total_ij != 0. + 0im
            Op_matrix[i,j] = total_ij
        end
    
    end
    return Op_matrix
end


function tuple_vector_equal(t::NTuple{N, Int}, v::Vector{Int}) where N
    @inbounds for i in 1:N
        if t[i] != v[i]
            return false
        end
    end
    return true
end

function calculate_matrix_elements_parallel(states::Vector{AbstractFockState}, Ops::MultipleFockOperator)
    n = length(states)
    Op_matrix = zeros(ComplexF64, n, n)

    Threads.@threads for idx in 1:(n*n)
        i = div(idx - 1, n) + 1
        j = mod(idx - 1, n) + 1
        bra = states[i]
        ket = states[j]
        total_ij = 0.0 + 0im
        tmp = MutableFockState(ket)
        
        for Op in Ops.terms
            reset2!(tmp, ket.occupations, ket.coefficient)
            apply!(Op, tmp)
            if tuple_vector_equal(bra.occupations, tmp.occupations)
                total_ij += bra.coefficient' * tmp.coefficient
            end
        end

        if total_ij != 0. + 0im
            Op_matrix[i, j] = total_ij
        end
    end
    return Op_matrix
end



function calculate_matrix_elements_naive(states::Vector{AbstractFockState}, Op::AbstractFockOperator)
    Op_matrix = zeros(ComplexF64, length(states), length(states))
    for (i,bra) in enumerate(states), (j,ket) in enumerate(states)  
        Op_matrix[i,j] = bra * (Op * ket)   
    end
    return Op_matrix
end

function sparseness(M::AbstractMatrix)
    s = 0
    t = 0
    for e in M 
        t+=1
        if e !=0.
            s+=1
        end
    end
    return s/t
end

end;
