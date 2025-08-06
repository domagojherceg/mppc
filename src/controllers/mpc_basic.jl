"""
    mpc_basic.jl
    
MPC controller with single control trajectory shared across all scenarios and single state trajectory.
"""

"""
    buildSMPC_RK4_basic(
        scf, scenV; 
        Ts = 1, 
        numWeeks = 12, 
        umax = 0.75, 
        T_max = 0.002, 
        rc = false, 
        i_alpha = 0.1,
        xf = false,
        singleCntrlTraj = false)
Build basic MPC controller with single control trajectory and initial condition.
"""

function buildSMPC_RK4_basic( scf, scenV; 
    Ts = 1, 
    numWeeks = 12, 
    umax = 0.75, 
    T_max = 0.002, 
    i_alpha = 0.1,
    xf = false)

    MPC = Model(Ipopt.Optimizer);
    set_attribute(MPC, "print_level", 0)

    numDays = numWeeks*7;
    numScenario = length(scenV)
    cntrlIntervals = Int64(numWeeks / Ts);
    numCnrlScenatio = length(scenV);
    
    #println(T_max)
    ## Optimization variables 
    @variable(MPC, S[1:numDays+1] ≥ 0 )
    @variable(MPC, I[1:numDays+1] ≥ 0 )
    @variable(MPC, D[1:numDays+1] ≥ 0 )
    @variable(MPC, T[1:numDays+1] ≥ 0 )
    @variable(MPC, H[1:numDays+1] ≥ 0 )
    @variable(MPC, E[1:numDays+1] ≥ 0 )
    @variable(MPC, 0 ≤ u[1:numDays] ≤ umax )



    ## initial point parameter
    @variable(MPC, p[i = 1:6] in Parameter(i))

    ## objective value 
    @objective(MPC, Min,  sum(scf.(u)) )

    ## enforce initial condition
    @constraint(MPC, S[1] .== p[1])
    @constraint(MPC, I[1] .== p[2])
    @constraint(MPC, D[1] .== p[3])
    @constraint(MPC, T[1] .== p[4])
    @constraint(MPC, H[1] .== p[5])
    @constraint(MPC, E[1] .== p[6]) 

        ## define RK4 discretrization 
    for idx_step in 1:numDays
        k1 = fsidthe(S[idx_step], 
               I[idx_step], 
               D[idx_step], 
               T[idx_step], 
               H[idx_step],
               E[idx_step], 
               scenV,
               u[idx_step])
        #possibly wrong, why is everyghing 1 here
        k2 = fsidthe(S[idx_step] + k1[1] / 2, 
               I[idx_step] + k1[2] / 2, 
               D[idx_step] + k1[3] / 2, 
               T[idx_step] + k1[4] / 2, 
               H[idx_step] + k1[5] / 2, 
               E[idx_step] + k1[6] / 2, 
               scenV, 
               u[idx_step])

        k3 = fsidthe(S[idx_step] + k2[1] / 2, 
               I[idx_step] + k2[2] / 2, 
               D[idx_step] + k2[3] / 2, 
               T[idx_step] + k2[4] / 2, 
               H[idx_step] + k2[5] / 2, 
               E[idx_step] + k2[6] / 2, 
               scenV, 
               u[idx_step])

        k4 = fsidthe(S[idx_step] + k3[1], 
               I[idx_step] + k3[2], 
               D[idx_step] + k3[3], 
               T[idx_step] + k3[4], 
               H[idx_step] + k3[5], 
               E[idx_step] + k3[6], 
               scenV,
               u[idx_step])

        @constraint(MPC, S[idx_step+1] == S[idx_step] + (1*k1[1] + 2*k2[1] + 2*k3[1] + 1*k4[1]) / 6.0)
        @constraint(MPC, I[idx_step+1] == I[idx_step] + (1*k1[2] + 2*k2[2] + 2*k3[2] + 1*k4[2]) / 6.0)
        @constraint(MPC, D[idx_step+1] == D[idx_step] + (1*k1[3] + 2*k2[3] + 2*k3[3] + 1*k4[3]) / 6.0)
        @constraint(MPC, T[idx_step+1] == T[idx_step] + (1*k1[4] + 2*k2[4] + 2*k3[4] + 1*k4[4]) / 6.0)
        @constraint(MPC, H[idx_step+1] == H[idx_step] + (1*k1[5] + 2*k2[5] + 2*k3[5] + 1*k4[5]) / 6.0)
        @constraint(MPC, E[idx_step+1] == E[idx_step] + (1*k1[6] + 2*k2[6] + 2*k3[6] + 1*k4[6]) / 6.0)
    end

    ## define the infection rate limit
    if rc == false
        @constraint(MPC, T .≤ T_max)
    end

    ## add constancy constraints on control input
    for idxInterval ∈ 1:cntrlIntervals
            @constraint(MPC, u[(idxInterval-1)*Ts*7 + 1] .== u[(idxInterval-1)*Ts*7 + 2 : (idxInterval-1)*Ts*7 + Ts*7]);
    end

    ## Xf constraints (terminal set)
    if xf == true
        for k = 1:numScenario
            tmax  = T_max;
            α = scenV[1]
            γ = scenV[2]
            λ = scenV[3]
            δ = scenV[4]
            σ = scenV[5]
            τ = scenV[6]
            np = (γ * λ) / (γ + λ)

            @constraint(MPC, S[numDays+1] ≤ (γ + np) / α / (1 -umax) )
            @constraint(MPC, I[numDays+1] ≤ (σ + τ) / δ  * (δ + λ) / γ * tmax)
            @constraint(MPC, D[numDays+1] ≤ (σ + τ) / δ * tmax)

        end
    end
    ## causality constraint
 
    return MPC
end