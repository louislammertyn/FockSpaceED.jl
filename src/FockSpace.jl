module FockSpace
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using VectorInterface
using KrylovKit
using SparseArrays
using Plots
using FFTW
using IterTools

include("./FockStates.jl")
include("./FockOps.jl")
include("./NormalOrder.jl")
include("./LatticeGeo.jl")
include("./TimeEv.jl")
include("./CommonOps.jl")
include("./utils.jl")
#####################################################################################################
#####################################################################################################


export AbstractFockSpace, U1FockSpace, UnrestrictedFockSpace,
       AbstractFockState, FockState, MultipleFockState, ZeroFockState
export fock_state, copy, cleanup_FS, checkU1, basisFS
export a_j, ad_j
export norm2FS
export all_states_U1, create_MFS, dot

export MutableFockState, to_fock_state, reset!, reset2!, norm2FS, cleanup_FS, mul_Mutable!
export a_j!, ad_j!

#####################################################################################################
#####################################################################################################


export AbstractFockOperator, FockOperator, MultipleFockOperator, ZeroFockOperator, identity_fockoperator
export calculate_matrix_elements,  calculate_matrix_elements_naive, calculate_matrix_elements_parallel, tuple_vector_equal
export sparseness
export cleanup_FO, dagger_FO
export apply!, apply
export rand_superpos, diagonalise_KR
export all_states_U1_O

#####################################################################################################
#####################################################################################################

export AbstractFockString, SameSiteString, MultiSiteString
export commute_first!, normal_order!, _normal_order!, normal_order, commutator

#####################################################################################################
#####################################################################################################

export vectorise_lattice, lattice_vectorisation_map, Lattice_NN, vector_to_lattice, Lattice, AbstractLattice
export Lattice_NN_h, vectorise_NN

#####################################################################################################
#####################################################################################################

export Time_Evolution, Time_Evolution_TD, schrodinger!, schrodinger_TD!

#####################################################################################################
#####################################################################################################

export MB_tensor, Entanglement_Entropy, density_onsite, one_body_ρ, density_flucs, momentum_density
export Bose_Hubbard_H, two_body_Op, four_body_Op, get_tensor_2body, get_tensor_4body, get_tensor_2body, delta, momentum_space_Op

#####################################################################################################
#####################################################################################################

export delta, nbody_geometry, fill_nbody_tensor, make_index
export periodic_neighbour, neighbour, helical
end
