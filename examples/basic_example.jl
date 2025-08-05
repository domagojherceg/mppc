"""
Example script for using the MPPC package with SIDTHE model.
This demonstrates basic usage of the package.
"""

using Pkg
# Activate the local package environment
Pkg.activate(".")

using MPPC
using Plots
using Printf
using Statistics

# Set parameters for simulation
numSteps = 2            # Number of MPC steps to simulate
predHorizonWeeks = 4     # Prediction horizon in weeks
Ts = 2                   # Sampling time in weeks (NPIs are kept constant for Ts weeks)
t_max = 0.002            # Maximum allowed infected proportion
use_multi_control = false  # Whether to use multiple control trajectories

# Get default parameters and generate scenarios
paramVal, probVal = get_default_parameters()
all_combinations, all_probabilities = genScenWithProbs(paramVal, probVal)
all_scenarios = eachrow(all_combinations)

# Choose cost function (squared cost)
fcost = squared_cost

# Set initial state
x_0 = get_default_initial_state()
println("Initial state: ", x_0)

# Build MPC controller based on configuration
if use_multi_control
    println("Using recourse MPC controller")
    mpc = buildSMPC_RK4_recourse(
        fcost, 
        all_scenarios, 
        paramVal,
        numWeeks = predHorizonWeeks,
        T_max = t_max,
        Ts = Ts,
        xf = true
    )
else
    println("Using robust MPC controller")
    mpc = buildSMPC_RK4_robust(
        fcost, 
        all_scenarios, 
        paramVal,
        numWeeks = predHorizonWeeks,
        T_max = t_max,
        Ts = Ts,
        xf = true
    )
end

# Choose scenario for simulation (using scenario 1)
scenario_index = 1
println("Running simulation with scenario: ", all_scenarios[scenario_index])

# Run MPC simulation
results = run_mpc_simulation(
    mpc,
    all_scenarios, 
    scenario_index, 
    numSteps, 
    x_0, 
    Ts = Ts
)

# Plot results
plots = plot_results(results, title = "SIDTHE MPC Simulation")
display(plots["combined"])

# Save plots if desired
savefig(plots["combined"], "sidthe_simulation.pdf")

# Print final results
println("\nSimulation completed.")
println("Final infected: ", results["trajectory_I"][end])
println("Mean control effort: ", mean(results["trajectory_U"]))
println("Mean objective value: ", mean(results["trajectory_J"]))
