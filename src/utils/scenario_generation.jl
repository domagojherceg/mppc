"""
    scenario_generation.jl
    
Utility functions for generating scenarios for stochastic MPC.
"""

"""
    genScen(a::Vector{T}, b::Vector{T}) where T

Generate all possible combinations of elements from vectors a and b.
Each element in the result is either a[i] or b[i].

Returns an array of unique result vectors.
"""
function genScen(a::Vector{T}, b::Vector{T}) where T
    result_vectors = []
    # c[i] is either a[i] or b[i]
    for combination in 0:1:(2^length(a) - 1)
        c = zeros(T, length(a))
        for i in 1:length(a)
            if (combination >> (i - 1)) & 1 != 0  # Test the i-th bit of combination
                c[i] = b[i]
            else
                c[i] = a[i]
            end
        end
        push!(result_vectors, c)
    end

    return unique(result_vectors)
end

"""
    genScenWithProbs(values::Matrix{T}, probMat::Matrix{T}) where T

Generate all possible combinations of parameter values with associated probabilities.

Input:
- values: Matrix where each row represents a parameter and each column represents a possible value
- probMat: Matrix of same dimensions as 'values' containing the probabilities for each value

Returns:
- all_combinations: Matrix where each row is a scenario (parameter combination)
- all_probabilities: Vector of scenario probabilities
"""
function genScenWithProbs(values::Matrix{T}, probMat::Matrix{T}) where T
    N, m = size(values)
    
    # Check that probabilities sum to 1 for each row
    for i in 1:N
        if abs(sum(probMat[i, :]) - 1.0) > 1e-10
            error("Probabilities in row $i do not sum to 1")
        end
    end
    
    # Generate all possible combinations
    num_combinations = m^N
    all_combinations = zeros(T, num_combinations, N)
    all_probabilities = zeros(T, num_combinations)
    
    for i in 0:(num_combinations - 1)
        combination = zeros(T, N)
        probability = 1.0
        for j in 1:N
            idx = div(i, m^(j - 1)) % m + 1
            combination[j] = values[j, idx]
            probability *= probMat[j, idx]
        end
        all_combinations[i + 1, :] = combination
        all_probabilities[i + 1] = probability
    end
    
    return all_combinations, all_probabilities
end
