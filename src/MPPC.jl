module MPPC

using DifferentialEquations, JuMP
using Ipopt, Distributions
using OrdinaryDiffEq, LinearAlgebra
using Printf, JLD2, Plots

# Models
include("models/sidthe_model.jl")

# Controllers
include("controllers/mpc_robust.jl")
include("controllers/mpc_recourse.jl")

# Utilities
include("utils/scenario_generation.jl")
include("utils/cost_functions.jl")
include("utils/simulation.jl")

# Export public API
export sidthe, 
       buildSMPC_RK4_robust, buildSMPC_RK4_recourse,
       genScenWithProbs, 
       squared_cost, linear_cost,
       run_mpc_simulation,
       plot_results,
       get_default_parameters, get_default_initial_state, get_imax

end # module
