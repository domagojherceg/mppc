"""
    simulation.jl
    
Utilities for running MPC simulations.
"""

"""
    run_mpc_simulation(
        mpc_model, 
        scenario_index, 
        num_steps, 
        initial_state; 
        Ts = 2,
        plot_results = false,
        debug_mode = false
    )

Run an MPC simulation for a given number of steps.

Parameters:
- mpc_model: MPC controller model
- scenario_index: Index of the scenario to use for simulation
- num_steps: Number of MPC steps to simulate
- initial_state: Initial state vector [S, I, D, T, H, E]
- Ts: Sampling time in weeks
- plot_results: Whether to plot results during simulation
- debug_mode: Whether to include additional debug information

Returns:
- Dictionary with trajectory data and MPC results
"""
function run_mpc_simulation(
    mpc_model, 
    all_scenarios, 
    scenario_index, 
    num_steps, 
    initial_state; 
    Ts = 2,
    plot_results = false,
    debug_mode = false
)
    # Initialize result containers
    trajectory_S = Vector{Float64}(undef, num_steps)
    trajectory_I = Vector{Float64}(undef, num_steps)
    trajectory_D = Vector{Float64}(undef, num_steps)
    trajectory_T = Vector{Float64}(undef, num_steps)
    trajectory_U = Vector{Float64}(undef, num_steps)
    trajectory_J = Vector{Float64}(undef, num_steps)
    
    # Detailed trajectories
    long_S = zeros(0)
    long_I = zeros(0)
    long_D = zeros(0)
    long_T = zeros(0)
    long_U = zeros(0)
    
    # MPC prediction outputs for all steps
    v_S_multi = Matrix{Float64}[]
    v_I_multi = Matrix{Float64}[]
    v_D_multi = Matrix{Float64}[]
    v_T_multi = Matrix{Float64}[]
    v_U_multi = Matrix{Float64}[]
    v_J_multi = Float64[]
    
    # Initialize state
    state = copy(initial_state)
    
    # Store initial state in detailed trajectories
    append!(long_S, state[1])
    append!(long_I, state[2])
    append!(long_D, state[3])
    append!(long_T, state[4])
    
    @printf " STEP | NPI   | SOLVER STATUS | TIME  |\n"
    @printf "------+-------+---------------+-------|\n"
    
    # Run MPC simulation
    for step = 1:num_steps
        # Store current state in trajectory
        trajectory_S[step] = state[1]
        trajectory_I[step] = state[2]
        trajectory_D[step] = state[3]
        trajectory_T[step] = state[4]
        
        # Update MPC model with current state
        set_parameter_value.(mpc_model[:p], state)
        
        # Solve MPC problem
        optimize!(mpc_model)
        
        # Store objective value
        trajectory_J[step] = objective_value(mpc_model)
        
        # Get optimal control input for the first time step
        if haskey(mpc_model, :u) && ndims(mpc_model[:u]) == 1
            # Single control trajectory
            npi = value(mpc_model[:u][1])
        elseif haskey(mpc_model, :u) && ndims(mpc_model[:u]) == 2
            # Multiple control trajectories (one per scenario)
            # Take first control from first scenario (they should be the same due to non-anticipativity)
            npi = value(mpc_model[:u][1, 1])
        else
            error("Unexpected control variable structure in MPC model")
        end
        
        trajectory_U[step] = npi
        
        # Store predicted trajectories
        push!(v_S_multi, value.(mpc_model[:S]))
        push!(v_I_multi, value.(mpc_model[:I]))
        push!(v_D_multi, value.(mpc_model[:D]))
        push!(v_T_multi, value.(mpc_model[:T]))
        
        if ndims(mpc_model[:u]) == 1
            # Single control trajectory
            push!(v_U_multi, reshape(value.(mpc_model[:u]), 1, :))
        else
            # Multiple control trajectories
            push!(v_U_multi, value.(mpc_model[:u]))
        end
        
        push!(v_J_multi, objective_value(mpc_model))
        
        # Simulate system for one sampling interval
        tspan = (0.0, Ts*7.0)
        prob = ODEProblem(sidthe, state, tspan, [all_scenarios[scenario_index]; npi])
        sol = solve(prob, saveat = 1, RK4())
        
        # Store detailed trajectories
        append!(long_S, sol[1, 2:end])
        append!(long_I, sol[2, 2:end])
        append!(long_D, sol[3, 2:end])
        append!(long_T, sol[4, 2:end])
        append!(long_U, repeat([npi], Ts*7))
        
        # Update state for next iteration
        state = sol[end]
        
        # Print progress
        @printf " %3d  | %.4f | %13s | %.2fs |\n" step npi termination_status(mpc_model) solve_time(mpc_model)
    end
    
    # Return results dictionary
    return Dict(
        "trajectory_S" => trajectory_S,
        "trajectory_I" => trajectory_I,
        "trajectory_D" => trajectory_D,
        "trajectory_T" => trajectory_T,
        "trajectory_U" => trajectory_U,
        "trajectory_J" => trajectory_J,
        "long_S" => long_S,
        "long_I" => long_I,
        "long_D" => long_D,
        "long_T" => long_T,
        "long_U" => long_U,
        "v_S_multi" => v_S_multi,
        "v_I_multi" => v_I_multi,
        "v_D_multi" => v_D_multi,
        "v_T_multi" => v_T_multi,
        "v_U_multi" => v_U_multi,
        "v_J_multi" => v_J_multi,
        "debug_info" => Dict(
            "scenario_used" => all_scenarios[scenario_index],
            "scenario_index" => scenario_index,
            "mpc_model" => mpc_model,
            "num_steps" => num_steps,
            "Ts" => Ts
        )
    )
end

"""
    plot_results(results; title="SIDTHE Model Simulation")

Plot the results of an MPC simulation.

Parameters:
- results: Dictionary with trajectory data from run_mpc_simulation
- title: Plot title

Returns:
- Dictionary of plots
"""
function plot_results(results; title="SIDTHE Model Simulation")
    p_I = plot(results["trajectory_I"], 
        linewidth=2,
        xlabel = "Time [weeks]",
        ylabel = "Infected",
        title = "Infected Population",
        label = "Infected",
        marker = :circle)
    
    p_U = plot(results["trajectory_U"], 
        linewidth=2,
        xlabel = "Time [weeks]",
        ylabel = "Control Input",
        title = "Control Input (NPI)",
        label = "u",
        marker = :circle)
    
    p_states = plot(1:length(results["trajectory_S"]), 
        [results["trajectory_S"] results["trajectory_I"] results["trajectory_D"] results["trajectory_T"]],
        linewidth=2,
        xlabel = "Time [weeks]",
        ylabel = "Population Fraction",
        title = "State Trajectories",
        label = ["Susceptible" "Infected" "Diagnosed" "Treated"],
        marker = :circle)
    
    p_long = plot(1:length(results["long_I"]),
        results["long_I"],
        linewidth=2,
        xlabel = "Time [days]",
        ylabel = "Population Fraction",
        title = "Daily Infected Trajectory",
        label = "Infected")
    
    # Combined plot
    p_combined = plot(p_states, p_U, p_I, p_long, 
        layout = (2, 2), 
        size = (1000, 800),
        plot_title = title)
    
    plots = Dict(
        "infected" => p_I,
        "control" => p_U,
        "states" => p_states,
        "long_infected" => p_long,
        "combined" => p_combined
    )
    
    return plots
end
