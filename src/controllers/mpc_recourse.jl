"""
    mpc_recourse.jl
    
MPC controller with multiple control trajectories (one for each scenario).
"""

"""
    buildSMPC_RK4_recourse(
        scf, scenV, scenP; 
        numWeeks = 12, 
        umax = 0.75, 
        T_max = 0.002, 
        Ts = 1,
        xf = false)

Build recourse stochastic MPC controller with scenario-specific control trajectories.

Parameters:
- scf: Cost function for control input
- scenV: Vector of scenarios (parameter sets)
- numWeeks: Prediction horizon in weeks
- umax: Maximum control input value
- T_max: Maximum allowed infected proportion
- Ts: Sampling time in weeks
- xf: Flag for terminal constraint

Returns:
- MPC model (JuMP model)
"""
function buildSMPC_RK4_recourse( scf, scenV, scenP; 
    Ts = 1, 
    numWeeks = 12, 
    umax = 0.75, 
    T_max = 0.002, 
    i_alpha = 0.1,
    xf = false,
    singleCntrlTraj = false,
    debug_mode = false)

    MPC = Model(Ipopt.Optimizer);
    set_attribute(MPC, "print_level", 0)

    numDays = numWeeks*7;
    numScenario = length(scenV)
    cntrlIntervals = Int64(numWeeks / Ts);
    numCnrlScenatio = length(scenV);
    
    #println(T_max)
    ## Optimization variables 
    @variable(MPC, S[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, I[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, D[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, T[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, H[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, E[1:numScenario,1:numDays+1] ≥ 0 )
    @variable(MPC, 0 ≤ u[1:numScenario,1:numDays] ≤ umax )

    ## initial point parameter
    @variable(MPC, p[i = 1:6] in Parameter(i))

    ## objective value 
    @objective(MPC, Min, sum(scenP[i] * sum(scf.(u[i,j]) for j in 1:numDays) for i in 1:numScenario) )

    ## enforce initial condition
    @constraint(MPC, S[:, 1] .== p[1])
    @constraint(MPC, I[:, 1] .== p[2])
    @constraint(MPC, D[:, 1] .== p[3])
    @constraint(MPC, T[:, 1] .== p[4])
    @constraint(MPC, H[:, 1] .== p[5])
    @constraint(MPC, E[:, 1] .== p[6]) 

        ## define RK4 discretrization 
    for k in 1:numScenario, idx_step in 1:numDays
        k1 = fsidthe(S[k, idx_step], 
               I[k, idx_step], 
               D[k, idx_step], 
               T[k, idx_step], 
               H[k, idx_step],
               E[k, idx_step], 
               scenV[k],
               u[k, idx_step])
        #possibly wrong, why is everyghing 1 here
        k2 = fsidthe(S[k, idx_step] + k1[1] / 2, 
               I[k, idx_step] + k1[2] / 2, 
               D[k, idx_step] + k1[3] / 2, 
               T[k, idx_step] + k1[4] / 2, 
               H[k, idx_step] + k1[5] / 2, 
               E[k, idx_step] + k1[6] / 2, 
               scenV[k], 
               u[k, idx_step])

        k3 = fsidthe(S[k, idx_step] + k2[1] / 2, 
               I[k, idx_step] + k2[2] / 2, 
               D[k, idx_step] + k2[3] / 2, 
               T[k, idx_step] + k2[4] / 2, 
               H[k, idx_step] + k2[5] / 2, 
               E[k, idx_step] + k2[6] / 2, 
               scenV[k], 
               u[k, idx_step])

        k4 = fsidthe(S[k, idx_step] + k3[1], 
               I[k, idx_step] + k3[2], 
               D[k, idx_step] + k3[3], 
               T[k, idx_step] + k3[4], 
               H[k, idx_step] + k3[5], 
               E[k, idx_step] + k3[6], 
               scenV[k],
               u[k, idx_step])

        @constraint(MPC, S[k, idx_step+1] == S[k, idx_step] + (1.0*k1[1] + 2.0*k2[1] + 2.0*k3[1] + 1.0*k4[1]) / 6.0)
        @constraint(MPC, I[k, idx_step+1] == I[k, idx_step] + (1.0*k1[2] + 2.0*k2[2] + 2.0*k3[2] + 1.0*k4[2]) / 6.0)
        @constraint(MPC, D[k, idx_step+1] == D[k, idx_step] + (1.0*k1[3] + 2.0*k2[3] + 2.0*k3[3] + 1.0*k4[3]) / 6.0)
        @constraint(MPC, T[k, idx_step+1] == T[k, idx_step] + (1.0*k1[4] + 2.0*k2[4] + 2.0*k3[4] + 1.0*k4[4]) / 6.0)
        @constraint(MPC, H[k, idx_step+1] == H[k, idx_step] + (1.0*k1[5] + 2.0*k2[5] + 2.0*k3[5] + 1.0*k4[5]) / 6.0)
        @constraint(MPC, E[k, idx_step+1] == E[k, idx_step] + (1.0*k1[6] + 2.0*k2[6] + 2.0*k3[6] + 1.0*k4[6]) / 6.0)
    end

    ## define the infection rate limit

    @constraint(MPC, T .≤ T_max)
  
    ## add constancy constraints on control input
    for k ∈ 1:numScenario, idxInterval ∈ 1:cntrlIntervals
            @constraint(MPC, u[k,(idxInterval-1)*Ts*7 + 1] .== u[k,(idxInterval-1)*Ts*7 + 2 : (idxInterval-1)*Ts*7 + Ts*7]);
    end

    ## Xf constraints (terminal set)
    if xf == true
        for k = 1:numScenario
            tmax  = T_max;
            α = scenV[k][1]
            γ = scenV[k][2]
            λ = scenV[k][3]
            δ = scenV[k][4]
            σ = scenV[k][5]
            τ = scenV[k][6]
            np = (γ * λ) / (γ + λ)

            @constraint(MPC, S[k, numDays+1] ≤ (γ + np) / α / (1 -umax) )
            @constraint(MPC, I[k, numDays+1] ≤ (σ + τ) / δ  * (δ + λ) / γ * tmax)
            @constraint(MPC, D[k, numDays+1] ≤ (σ + τ) / δ * tmax)

        end
    end
    ## causality constraint
    @constraint(MPC, u[1,1] .== u[2:end,1]);

    # Return MPC model with or without debug info
    if debug_mode
        return MPC, Dict(
            "S" => S,
            "I" => I,
            "D" => D,
            "T" => T,
            "H" => H,
            "E" => E,
            "u" => u,
            "numDays" => numDays,
            "numScenario" => numScenario
        )
    else
        return MPC
    end
end
