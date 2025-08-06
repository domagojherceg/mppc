"""
    cost_functions.jl
    
Various cost functions for MPC optimization.
"""

"""
    squared_cost(u)

Squared cost function for control input u.
"""
function squared_cost(u)
    return u^2
end

"""
    linear_cost(u)

Linear cost function for control input u.
"""
function linear_cost(u)
    return u
end
