using Statistics
using StatsBase
using Distributions
using LinearAlgebra
using Random
using Printf

# ==========================================
# 1. DBBN Model Definition
# ==========================================

"""
DiscreteDBBN Struct
Stores the Conditional Probability Table (CPT) as a transition matrix
and the representative values for each discrete state.
"""
struct DiscreteDBBN
    state_centers::Vector{Float64}      # Physical value for each state (e.g., wind speed m/s)
    transition_matrix::Matrix{Float64}  # P(State_t | State_{t-1})
    initial_dist::Vector{Float64}       # Initial state probability P(State_0)
end

# ==========================================
# 2. Helper: Generate "Mock" Historical Data
# ==========================================

function generate_mock_history_data(n_hours::Int)
    # Simulate Wind Speed (Weibull + Autocorrelation)
    wind_speed = zeros(n_hours)
    wind_speed[1] = 5.0
    for t in 2:n_hours
        # AR(1) process
        noise = rand(Normal(0, 1.5))
        wind_speed[t] = 0.9 * wind_speed[t-1] + 0.1 * 7.0 + noise
        wind_speed[t] = clamp(wind_speed[t], 0.0, 25.0)
    end
    
    # Simulate Solar Irradiance (Daily Cycle + Cloud noise)
    solar_irrad = zeros(n_hours)
    for t in 1:n_hours
        hour = (t - 1) % 24
        if 6 <= hour <= 18
            # Max theoretical * random cloud factor
            max_solar = 1000.0 * sin((hour - 6) * π / 12)
            cloud_factor = rand(Beta(2, 1)) # Skewed towards sunny
            solar_irrad[t] = max_solar * cloud_factor
        else
            solar_irrad[t] = 0.0
        end
    end
    
    return wind_speed, solar_irrad
end

# ==========================================
# 3. DBBN Core: Training (MLE)
# ==========================================

"""
Discretize continuous data into state indices
"""
function discretize_data(data::Vector{Float64}, n_states::Int)
    min_val, max_val = minimum(data), maximum(data)
    
    # Use range to avoid floating point issues with step addition
    # Create n_states bins
    edges = range(min_val, max_val, length=n_states+1)
    
    state_centers = [(edges[i] + edges[i+1])/2 for i in 1:n_states]
    
    discrete_seq = Int[]
    sizehint!(discrete_seq, length(data))
    
    for x in data
        # Find bin index
        idx = searchsortedlast(edges, x)
        # Handle the case where x == min_val (idx returns 1, correct)
        # Handle the case where x == max_val (idx returns n_states+1, needs clamping)
        idx = clamp(idx, 1, n_states)
        push!(discrete_seq, idx)
    end
    
    return discrete_seq, state_centers
end

"""
Train DBBN: Calculate transition matrix using Maximum Likelihood Estimation (MLE)
Corresponds to Section II.B in the paper
"""
function train_dbbn(data::Vector{Float64}, n_states::Int=10)
    # 1. Discretize
    seq, centers = discretize_data(data, n_states)
    
    # 2. Initialize counts
    counts = zeros(Int, n_states, n_states) # counts[i, j]: transition i -> j
    state_counts = zeros(Int, n_states)     # Total occurrences of each state
    
    # 3. Count transitions
    for t in 2:length(seq)
        prev_s = seq[t-1]
        curr_s = seq[t]
        counts[prev_s, curr_s] += 1
        state_counts[prev_s] += 1 # Count occurrence for initial dist approximation
    end
    # Add last state to total counts
    state_counts[seq[end]] += 1
    
    # 4. Normalize to probabilities
    trans_mat = zeros(Float64, n_states, n_states)
    for i in 1:n_states
        row_sum = sum(counts[i, :])
        if row_sum > 0
            trans_mat[i, :] = counts[i, :] ./ row_sum
        else
            # Absorbing state or unseen transition: stay in same state
            trans_mat[i, i] = 1.0 
        end
    end
    
    # Estimate initial distribution using stationary frequency (histogram)
    # This is more robust than just taking the first data point
    total_obs = sum(state_counts)
    init_dist = state_counts ./ total_obs
    
    # Fallback if something went wrong (shouldn't happen with valid data)
    if sum(init_dist) == 0
        init_dist = fill(1.0/n_states, n_states)
    end
    
    return DiscreteDBBN(centers, trans_mat, init_dist)
end

# ==========================================
# 4. DBBN Core: Sampling/Prediction
# ==========================================

"""
Generate new time series using the trained DBBN
"""
function sample_dbbn(model::DiscreteDBBN, n_steps::Int, start_val::Float64=NaN)
    n_states = length(model.state_centers)
    trajectory = zeros(Float64, n_steps)
    
    # Determine initial state
    current_state = 1
    if isnan(start_val)
        # Sample from initial distribution
        current_state = rand(Categorical(model.initial_dist))
    else
        # Find state closest to start_val
        current_state = argmin(abs.(model.state_centers .- start_val))
    end
    
    trajectory[1] = model.state_centers[current_state]
    
    # Step-by-step sampling
    for t in 2:n_steps
        # Get transition probabilities for current state
        probs = model.transition_matrix[current_state, :]
        
        # Ensure probabilities sum to 1 (fix floating point errors)
        probs = probs ./ sum(probs)
        
        # Sample next state
        next_state = rand(Categorical(probs))
        
        # Record value (add small noise to smooth discretization steps)
        val = model.state_centers[next_state]
        # Optional: Add small smoothing noise
        # val += (rand() - 0.5) * (model.state_centers[2] - model.state_centers[1]) * 0.5
        
        trajectory[t] = val
        current_state = next_state
    end
    
    return trajectory
end

# ==========================================
# 5. Integration: Physical Conversion
# ==========================================

function wind_power_curve(speed_ms)
    v_in = 3.0
    v_rated = 12.0
    v_out = 25.0
    p_rated = 450.0 # MW
    
    if speed_ms < v_in || speed_ms > v_out
        return 0.0
    elseif speed_ms >= v_rated
        return p_rated
    else
        return p_rated * ((speed_ms - v_in) / (v_rated - v_in))^3
    end
end

# ==========================================
# 6. Demo
# ==========================================

function run_dbbn_demo()
    println("1. Generating Mock Historical Data (Training Set)...")
    hist_wind, hist_solar = generate_mock_history_data(8760)
    
    println("2. Training DBBN Models (MLE)...")
    wind_dbbn = train_dbbn(hist_wind, 15)
    
    # Note: Simple Markov Chain is poor for Solar due to diurnal cycles.
    # In a real implementation, we would model the "Clearness Index" or use a non-homogeneous model.
    solar_dbbn = train_dbbn(hist_solar, 10)
    
    println("3. Generating New Scenarios using DBBN...")
    T_horizon = 24
    
    # Sample Wind
    sim_wind_speed = sample_dbbn(wind_dbbn, T_horizon, hist_wind[end])
    sim_wind_power = wind_power_curve.(sim_wind_speed)
    
    # Sample Solar
    sim_solar_power = sample_dbbn(solar_dbbn, T_horizon, hist_solar[end])
    
    # Post-process Solar: Force night to zero
    for t in 1:T_horizon
        hour_of_day = (t - 1) % 24
        if hour_of_day < 6 || hour_of_day > 18
            sim_solar_power[t] = 0.0
        end
    end

    println("\n--- DBBN Simulation Result (First 10 hours) ---")
    println("Hour | Wind Speed (m/s) | Wind Power (MW) | Solar Power (MW)")
    for t in 1:10
        @printf("%4d | %16.2f | %15.2f | %16.2f\n", 
            t, sim_wind_speed[t], sim_wind_power[t], sim_solar_power[t])
    end
    
    return sim_wind_power, sim_solar_power
end

"""
    twist_dbbn(model::DiscreteDBBN, twisting_factor::Float64)

Applies exponential twisting to the transition matrix to favor specific states.
If twisting_factor > 0, it favors LOWER value states (stress testing for low wind/solar).
If twisting_factor < 0, it favors HIGHER value states.

Returns:
    1. A new DiscreteDBBN with the biased transition matrix.
    2. The original model (for likelihood ratio calculation).
"""
function twist_dbbn(model::DiscreteDBBN, twisting_factor::Float64)
    n_states = length(model.state_centers)
    new_trans = zeros(Float64, n_states, n_states)
    
    # We assume state indices correlate with magnitude (1 is low, N is high)
    # We want to bias towards state 1.
    # Weight for state j: w_j = exp(-twisting_factor * j)
    
    for i in 1:n_states
        denom = 0.0
        # Calculate normalization factor for this row
        for j in 1:n_states
            # Original prob * bias weight
            # Using normalized state index (0 to 1) for stability
            bias = exp(-twisting_factor * (j / n_states))
            new_trans[i, j] = model.transition_matrix[i, j] * bias
            denom += new_trans[i, j]
        end
        
        # Normalize row
        if denom > 0
            new_trans[i, :] ./= denom
        else
            new_trans[i, i] = 1.0
        end
    end
    
    # Also twist initial distribution
    new_init = zeros(Float64, n_states)
    denom_init = 0.0
    for j in 1:n_states
        bias = exp(-twisting_factor * (j / n_states))
        new_init[j] = model.initial_dist[j] * bias
        denom_init += new_init[j]
    end
    new_init ./= denom_init
    
    return DiscreteDBBN(model.state_centers, new_trans, new_init)
end

"""
    sample_dbbn_with_likelihood(original::DiscreteDBBN, biased::DiscreteDBBN, n_steps::Int)

Samples a trajectory from the BIASED model, but calculates the Likelihood Ratio (Weight).
Weight = P_original(Trajectory) / P_biased(Trajectory)
"""
function sample_dbbn_with_likelihood(original::DiscreteDBBN, biased::DiscreteDBBN, n_steps::Int)
    trajectory = zeros(Float64, n_steps)
    log_likelihood_ratio = 0.0
    
    # 1. Sample Initial State from BIASED distribution
    current_state = rand(Categorical(biased.initial_dist))
    trajectory[1] = biased.state_centers[current_state]
    
    # Update Weight: P_orig(s0) / P_biased(s0)
    p_orig = original.initial_dist[current_state]
    p_bias = biased.initial_dist[current_state]
    
    # Avoid log(0)
    if p_bias > 0 && p_orig > 0
        log_likelihood_ratio += log(p_orig) - log(p_bias)
    end
    
    # 2. Step-by-step sampling
    for t in 2:n_steps
        # Get transition probabilities from BIASED model
        probs_bias = biased.transition_matrix[current_state, :]
        
        # Sample next state using BIASED probabilities
        next_state = rand(Categorical(probs_bias))
        
        # Calculate probabilities of this specific transition in both models
        p_trans_orig = original.transition_matrix[current_state, next_state]
        p_trans_bias = biased.transition_matrix[current_state, next_state]
        
        # Update Weight
        if p_trans_bias > 0 && p_trans_orig > 0
            log_likelihood_ratio += log(p_trans_orig) - log(p_trans_bias)
        end
        
        trajectory[t] = biased.state_centers[next_state]
        current_state = next_state
    end
    
    # Return trajectory and the actual Likelihood Ratio (exp(log_LR))
    return trajectory, exp(log_likelihood_ratio)
end

# Run Demo
sim_wind, sim_solar = run_dbbn_demo()

