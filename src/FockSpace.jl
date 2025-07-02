module FockSpace
using LinearAlgebra
using DifferentialEquations

include("./FockStates.jl")
include("./FockOps.jl")
include("./LatticeGeo.jl")
include("./TimeEv.jl")
include("./CommonOps.jl")
#####################################################################################################
#####################################################################################################


export AbstractFockSpace, U1FockSpace, UnrestrictedFockSpace,
       AbstractFockState, FockState, MultipleFockState, ZeroFockState
export fock_state, copy, cleanup_FS, checkU1
export a_j, ad_j
export norm2FS
export all_states_U1, create_MFS

export  MutableFockState, to_fock_state, reset!, reset2!, norm2FS, cleanup_FS, mul_Mutable!
export a_j!, ad_j!

#####################################################################################################
#####################################################################################################


export AbstractFockOperator, FockOperator, MultipleFockOperator, ZeroFockOperator
export calculate_matrix_elements,  calculate_matrix_elements_naive, calculate_matrix_elements_parallel, tuple_vector_equal
export sparseness
export cleanup_FO, dagger_FO
export apply!

#####################################################################################################
#####################################################################################################

export vectorise_lattice, lattice_vectorisation_map, Lattice_NN

#####################################################################################################
#####################################################################################################

export Time_Evolution, Time_Evolution_TD, schrodinger!, schrodinger_TD!

#####################################################################################################
#####################################################################################################

export MB_tensor, Entanglement_Entropy, density_onsite, one_body_ρ, density_flucs, momentum_density

end
