# MPPC

A Julia package for Model Predictive Control (MPC) with epidemiological models, specifically the SIDTHE (Susceptible, Infected, Diagnosed, Treated, Healed, Extinct) model.

## Features

- Stochastic Model Predictive Control (SMPC) for epidemic control
- Single and multiple control trajectory approaches
- Visualization and simulation utilities

## Directory Structure

```
MPPC/
├── examples/
│   ├── basic_example.jl
│   └── compare_control_strategies.jl
├── src/
│   ├── controllers/
│   │   ├── mpc_basic.jl
│   │   ├── mpc_recourse.jl
│   │   └── mpc_robust.jl
│   ├── models/
│   │   └── sidthe_model.jl
│   ├── utils/
│   │   ├── cost_functions.jl
│   │   ├── scenario_generation.jl
│   │   └── simulation.jl
│   └── MPPC.jl
├── test/
└── Project.toml
```

## Installation

This is currently a local Julia package. To use it:

1. Clone the repository
2. Navigate to the project directory
3. Start Julia in the project directory
4. Activate the environment:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Usage

Basic example:

```julia
using MPPC

# Generate scenarios
paramVal, probVal = get_default_parameters()
all_combinations, all_probabilities = genScenWithProbs(paramVal, probVal)
all_scenarios = eachrow(all_combinations)

# Create MPC controller
mpc = buildSMPC_RK4_robust(
    squared_cost,
    all_scenarios,
    paramVal,
    numWeeks = 4,
    T_max = 0.002,
    Ts = 2,
    xf = true
)

# Run simulation
results = run_mpc_simulation(
    mpc,
    all_scenarios,
    1,  # scenario index
    10, # number of steps
    get_default_initial_state()
)

# Plot results
plots = plot_results(results)
```

## Examples

Check the `examples/` directory for more detailed examples:

- `basic_example.jl`: Simple demonstration of the package
- `compare_control_strategies.jl`: Comparison of single vs. multi-control approaches

## License

MIT
