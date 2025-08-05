"""
    sidthe_model.jl
    
Contains the SIDTHE (Susceptible, Infected, Diagnosed, Treated, Healed, Extinct) model
for pandemic modeling and control.
"""

"""
    sidthe(xcurr, p, t)

The SIDTHE model differential equation function.
- xcurr: State vector [S, I, D, T, H, E]
- p: Parameter vector [α, γ, λ, δ, σ, τ, u]
- t: Time (not used, included for compatibility with ODE solvers)

Returns the derivatives [dS/dt, dI/dt, dD/dt, dT/dt, dH/dt, dE/dt]
"""
function sidthe(xcurr, p, t)
    x = xcurr[1:end]
    α = p[1]    # Infection rate
    γ = p[2]    # Diagnosis rate
    λ = p[3]    # Death rate (not diagnosed)
    δ = p[4]    # Treatment rate
    σ = p[5]    # Recovery rate (from treatment)
    τ = p[6]    # Death rate (from treatment)
    
    # Non-pharmaceutical intervention
    u = length(p) > 6 ? p[end] : 0.0
    
    np = (γ*λ) / (γ + λ)
    
    return [
        -α*(1 - u)*x[1]*x[2],                  # dS/dt
        α*(1 - u)*x[1]*x[2] - (γ + np)*x[2],   # dI/dt 
        γ*x[2] - (δ+λ)*x[3],                   # dD/dt 
        δ*x[3] - (σ+τ)*x[4],                   # dT/dt
        σ*x[4] + λ*x[3]+ np*x[2],              # dH/dt
        τ*x[4]                                 # dE/dt        
    ]
end

"""
    fsidthe(s, i, d, t, h, e, p, u)

The SIDTHE model function for use in JuMP models.
- s, i, d, t, h, e: Individual state variables
- p: Parameter vector [α, γ, λ, δ, σ, τ]
- u: Control input (non-pharmaceutical intervention)

Returns the derivatives [dS/dt, dI/dt, dD/dt, dT/dt, dH/dt, dE/dt]
"""
function fsidthe(s, i, d, t, h, e, p, u)
    α = p[1]    # Infection rate
    γ = p[2]    # Diagnosis rate
    λ = p[3]    # Death rate (not diagnosed)
    δ = p[4]    # Treatment rate
    σ = p[5]    # Recovery rate
    τ = p[6]    # Death rate 
    
    np = (γ*λ) / (γ + λ)

    return [
        -α*(1 - u)*s*i,                 # dS/dt
        α*(1 - u)*s*i - (γ + np)*i,     # dI/dt 
        γ*i - (δ+λ)*d,                  # dD/dt 
        δ*d - (σ+τ)*t,                  # dT/dt
        σ*t + λ*d + np*i,               # dH/dt
        τ*t                             # dE/dt        
    ]
end

"""
    get_default_parameters()

Returns default parameter values for the SIDTHE model.
- paramMid: Middle parameter values
- paramLo: Lower parameter values
- paramHi: Higher parameter values
- probVal: Probabilities for each parameter set
"""
function get_default_parameters()
    # Scaling factors for low and high parameter values
    cl = 0.95
    cu = 1.05
    
    # Mid-range parameter values [α, γ, λ, δ, σ, τ]
    paramMid = [0.35, 0.1, 0.09, 0.002, 0.015, 0.01]
    
    # Low and high parameter values
    paramLo = cl .* paramMid
    paramHi = cu .* paramMid
    
    # Parameter value matrix (columns are scenarios)
    paramVal = [paramLo paramMid paramHi]
    
    # Probabilities for each parameter set
    probLo = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    probMid = [0.90, 0.90, 0.90, 0.90, 0.90, 0.90]
    probHi = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    
    probVal = [probLo probMid probHi]
    
    return paramVal, probVal
end

"""
    get_default_initial_state()

Returns the default initial state for SIDTHE model [S, I, D, T, H, E]
"""
function get_default_initial_state()
    return [0.99, 0.008, 0.0019, 0.00010, 0.000, 0.0]
end

"""
    get_imax()

Returns default maximum allowed infected proportion
"""
function get_imax()
    return 0.002
end
