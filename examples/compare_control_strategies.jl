"""
Example script for comparing single vs multi-control MPC strategies
for the SIDTHE model.
"""

using Pkg
# Activate the local package environment
Pkg.activate(".")

using MPPC
using Plots
using Printf
using Statistics
# Set parameters for simulation
numSteps = 3             # Number of MPC steps to simulate
predHorizonWeeks = 4     # Prediction horizon in weeks
Ts = 2                   # Sampling time in weeks
Tmax = 0.002             # Maximum allowed infected proportion
risk_alpha = 1           # Risk parameter

# Get default parameters and generate scenarios
paramVal, probVal = get_default_parameters()
all_combinations, all_probabilities = genScenWithProbs(paramVal, probVal)
all_scenarios = collect(eachrow(all_combinations))

# Choose cost function (squared cost)
fcost = squared_cost

# Set initial state
x_0 = get_default_initial_state()
println("Initial state: ", x_0)

# Build robust MPC controller
println("\n== Building robust MPC controller ==")
mpc_robust = buildSMPC_RK4_robust(
    fcost, 
    all_scenarios, 
    all_probabilities,
    numWeeks = predHorizonWeeks,
    T_max = Tmax,
    Ts = Ts,
    xf = true
)

# Build recourse MPC controller
println("\n== Building recourse MPC controller ==")
mpc_recourse = buildSMPC_RK4_recourse(
    fcost, 
    all_scenarios, 
    all_probabilities,
    numWeeks = predHorizonWeeks,
    T_max = Tmax,
    Ts = Ts,
    xf = true
)

# Choose scenario for simulation (using scenario 1)
scenario_index = 1
println("Running simulation with scenario: ", all_scenarios[scenario_index])
println("Initial state: ", x_0)
# Run MPC simulation with robust control
println("\n== Running simulation with robust MPC ==")
results_robust = run_mpc_simulation(
    mpc_robust,
    all_scenarios, 
    scenario_index, 
    numSteps, 
    x_0, 
    Ts = Ts
)

# Run MPC simulation with recourse control
println("\n== Running simulation with recourse MPC ==")
results_recourse = run_mpc_simulation(
    mpc_recourse,
    all_scenarios, 
    scenario_index, 
    numSteps, 
    x_0, 
    Ts = Ts
)
# Plot results for robust control
plots_robust = plot_results(results_robust, title = "SIDTHE MPC - Robust Control")

# Plot results for recourse control
plots_recourse = plot_results(results_recourse, title = "SIDTHE MPC - Recourse Control")

# Create comparison plots
p_I_compare = plot(
    1:numSteps,
    [results_robust["trajectory_I"] results_recourse["trajectory_I"]],
    linewidth=2,
    xlabel = "Time [weeks]",
    ylabel = "Infected",
    title = "Infected Population Comparison",
    label = ["Robust Control" "Recourse Control"],
    marker = :circle
)

p_U_compare = plot(
    1:numSteps,
    [results_robust["trajectory_U"] results_recourse["trajectory_U"]],
    linewidth=2,
    xlabel = "Time [weeks]",
    ylabel = "Control Input",
    title = "Control Input (NPI) Comparison",
    label = ["Robust Control" "Recourse Control"],
    marker = :circle
)

p_J_compare = plot(
    1:numSteps,
    [results_robust["trajectory_J"] results_recourse["trajectory_J"]],
    linewidth=2,
    xlabel = "Time [weeks]",
    ylabel = "Objective Value",
    title = "Objective Value Comparison",
    label = ["Robust Control" "Recourse Control"],
    marker = :circle
)

# Combined comparison plot
p_combined_compare = plot(
    p_I_compare, p_U_compare, p_J_compare,
    layout = (3, 1), 
    size = (800, 900),
    plot_title = "Single vs Multi Control Comparison"
)

display(p_combined_compare)

# Save comparison plot
savefig(p_combined_compare, "control_comparison.pdf")

# Print final results comparison
println("\n==== Simulation Results Comparison ====")
println("\nRobust Control:")
println("  Final infected: ", results_robust["trajectory_I"][end])
println("  Mean control effort: ", mean(results_robust["trajectory_U"]))
println("  Mean objective value: ", mean(results_robust["trajectory_J"]))

println("\nRecourse Control:")
println("  Final infected: ", results_recourse["trajectory_I"][end])
println("  Mean control effort: ", mean(results_recourse["trajectory_U"]))
println("  Mean objective value: ", mean(results_recourse["trajectory_J"]))

# Calculate improvements
i_improvement = (results_robust["trajectory_I"][end] - results_recourse["trajectory_I"][end]) / results_robust["trajectory_I"][end] * 100
u_improvement = (mean(results_robust["trajectory_U"]) - mean(results_recourse["trajectory_U"])) / mean(results_robust["trajectory_U"]) * 100
j_improvement = (mean(results_robust["trajectory_J"]) - mean(results_recourse["trajectory_J"])) / mean(results_robust["trajectory_J"]) * 100

println("\nImprovements (Recourse vs Robust):")
println("  Infected reduction: ", round(i_improvement, digits=2), "%")
println("  Control effort reduction: ", round(u_improvement, digits=2), "%")
println("  Objective value improvement: ", round(j_improvement, digits=2), "%")


