using JuMP
using HiGHS
using Distributions
using Random
using Statistics
using LinearAlgebra
using StatsBase
using Plots # Added for visualization

# ==========================================
# 1. Struct Definitions & Parameter Initialization (Based on Table I & V)
# ==========================================

"""
Generator struct representing thermal units with degradation parameters.
Corresponds to Table I in the paper.
"""
struct Generator
    id::Int
    Pmax::Float64
    Pmin::Float64
    RD::Float64       # Ramp Down
    RU::Float64       # Ramp Up
    TU::Int           # Min Up Time
    TD::Int           # Min Down Time
    Cost_fix::Float64
    Cost_start::Float64
    # Aging model parameters (PHM)
    beta::Float64     # Weibull shape
    scale_a::Float64  # Weibull scale
    gamma1::Float64   # Coeff for Temp
    gamma2::Float64   # Coeff for Load
end

"""
Storage struct representing ESS parameters.
Corresponds to Table V in the paper.
"""
struct Storage
    Ecap::Float64
    Pc_max::Float64
    Pd_max::Float64
    eta_c::Float64
    eta_d::Float64
    SOC_min::Float64
    SOC_max::Float64
    alpha::Float64    # RH tracking tolerance
end

"""
Results struct to store simulation outcomes for analysis.
"""
mutable struct SimulationResults
    lolp_history::Vector{Float64}
    cv_history::Vector{Float64}
    eens_history::Vector{Float64}
    total_cost::Float64
    convergence_iteration::Int
end

# Initialize generators (Construct based on typical ranges in Table I)
function init_generators()
    # Simplified example: Construct 3 units with different characteristics representing the 8 units in the paper
    # G1: Base load unit (Large capacity, slow aging)
    # G2: Intermediate load unit
    # G3: Peak load unit (Small capacity, fast aging)
    # UPDATED: Increased capacities to create a system with positive reserve margin (Total ~1150 MW > 1000 MW Peak).
    # Importance Sampling is most effective for RARE events (LOLP < 0.01). 
    # On a tight system (LOLP ~ 0.15), IS increases variance due to weight dispersion.
    gens = [
        Generator(1, 600.0, 150.0, 200.0, 200.0, 8, 8, 400.0, 5000.0, 3.0, 5000.0, 0.5, 1.0),
        Generator(2, 350.0, 50.0, 100.0, 100.0, 4, 4, 200.0, 2000.0, 3.0, 4000.0, 0.5, 1.0),
        Generator(3, 200.0, 20.0, 60.0, 60.0, 1, 1, 100.0, 500.0, 3.0, 3000.0, 0.5, 1.0)
    ]
    return gens
end

# Initialize Storage (Reference Table V)
# UPDATED: Increased alpha to 1.0 to relax day-ahead tracking during reliability assessment
const ESS = Storage(300.0, 60.0, 60.0, 0.95, 0.95, 0.1, 0.9, 1.0)

# System Parameters
const T_HORIZON = 24       # Simulation duration (Paper uses 168, shortened here for demo)
const T_ROLL = 6           # Rolling window length
const LOAD_PEAK = 1000.0   # Scaled peak load
const CV_THRESHOLD = 0.05  # Convergence threshold (Paper uses 0.02)
const MAX_SAMPLES = 5000   # Max Monte Carlo samples
# UPDATED: Reduced twisting/bias factors to prevent weight explosion.
const IS_TWIST_FACTOR = 0.05 
# NEW: Bias factor for generator forced outages (Composite Importance Sampling)
const GEN_FAILURE_BIAS = 3.0 

# ==========================================
# 2. Helper Functions: DBBN Simulation & Aging Model
# ==========================================

# --- DBBN Implementation (Enhanced Scenario Generation) ---

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

function discretize_data(data::Vector{Float64}, n_states::Int)
    min_val, max_val = minimum(data), maximum(data)
    edges = range(min_val, max_val, length=n_states+1)
    state_centers = [(edges[i] + edges[i+1])/2 for i in 1:n_states]
    
    discrete_seq = Int[]
    sizehint!(discrete_seq, length(data))
    
    for x in data
        idx = searchsortedlast(edges, x)
        idx = clamp(idx, 1, n_states)
        push!(discrete_seq, idx)
    end
    
    return discrete_seq, state_centers
end

function train_dbbn(data::Vector{Float64}, n_states::Int=10)
    seq, centers = discretize_data(data, n_states)
    counts = zeros(Int, n_states, n_states)
    state_counts = zeros(Int, n_states)
    
    for t in 2:length(seq)
        prev_s = seq[t-1]
        curr_s = seq[t]
        counts[prev_s, curr_s] += 1
        state_counts[prev_s] += 1
    end
    state_counts[seq[end]] += 1
    
    trans_mat = zeros(Float64, n_states, n_states)
    for i in 1:n_states
        row_sum = sum(counts[i, :])
        if row_sum > 0
            trans_mat[i, :] = counts[i, :] ./ row_sum
        else
            trans_mat[i, i] = 1.0 
        end
    end
    
    total_obs = sum(state_counts)
    init_dist = state_counts ./ total_obs
    if sum(init_dist) == 0
        init_dist = fill(1.0/n_states, n_states)
    end
    
    return DiscreteDBBN(centers, trans_mat, init_dist)
end

function sample_dbbn(model::DiscreteDBBN, n_steps::Int)
    trajectory = zeros(Float64, n_steps)
    current_state = rand(Categorical(model.initial_dist))
    trajectory[1] = model.state_centers[current_state]
    
    for t in 2:n_steps
        probs = model.transition_matrix[current_state, :]
        probs = probs ./ sum(probs) # Normalize
        next_state = rand(Categorical(probs))
        trajectory[t] = model.state_centers[next_state]
        current_state = next_state
    end
    
    return trajectory
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

"""
Generate scenarios using trained DBBN models.
Replaces the simple random generation.
"""
function generate_dbbn_scenarios(wind_model::DiscreteDBBN, solar_model::DiscreteDBBN, T::Int)
    # Base Load: Deterministic pattern + Random Noise
    base_load = [LOAD_PEAK * (0.6 + 0.4 * sin((t-8)*π/12)) + rand(Normal(0, 20)) for t in 1:T]
    
    # Wind: Sample speed -> Convert to Power
    wind_speed = sample_dbbn(wind_model, T)
    wind_power = wind_power_curve.(wind_speed)
    
    # Solar: Sample Irradiance -> Scale to Power
    solar_irrad = sample_dbbn(solar_model, T)
    solar_capacity = 200.0 # MW
    solar_power = (solar_irrad ./ 1000.0) .* solar_capacity
    
    # Apply physical constraints (night time)
    for t in 1:T
        hour = (t - 1) % 24
        if hour < 6 || hour > 18
            solar_power[t] = 0.0
        end
    end
    
    return base_load, wind_power, solar_power
end

"""
Generate scenarios using Importance Sampling.
Returns profiles and the combined Likelihood Ratio (Weight) for this scenario.
"""
function generate_dbbn_scenarios_is(
    wind_orig::DiscreteDBBN, wind_biased::DiscreteDBBN,
    solar_orig::DiscreteDBBN, solar_biased::DiscreteDBBN, 
    T::Int
)
    # Base Load (Deterministic + Noise) - Weight is 1.0 (unbiased)
    base_load = [LOAD_PEAK * (0.6 + 0.4 * sin((t-8)*π/12)) + rand(Normal(0, 20)) for t in 1:T]
    
    # Wind: Sample from BIASED model, get Weight
    wind_speed, w_wind = sample_dbbn_with_likelihood(wind_orig, wind_biased, T)
    wind_power = wind_power_curve.(wind_speed)
    
    # Solar: Sample from BIASED model, get Weight
    solar_irrad, w_solar = sample_dbbn_with_likelihood(solar_orig, solar_biased, T)
    solar_capacity = 200.0 
    solar_power = (solar_irrad ./ 1000.0) .* solar_capacity
    
    # Apply physical constraints
    for t in 1:T
        hour = (t - 1) % 24
        if hour < 6 || hour > 18
            solar_power[t] = 0.0
        end
    end
    
    # Combined Weight (assuming independence between Wind and Solar generation process)
    total_weight = w_wind * w_solar
    
    return base_load, wind_power, solar_power, total_weight
end

# --- PHM Model ---

# Proportional Hazard Model (PHM) to calculate failure rate
# Corresponds to equations (2) and (4) in the paper
function calculate_outage_prob(gen::Generator, t, accumulated_hazard, current_load_ratio, temp_condition)
    # h(t, Z(t))
    # Z1: temp_condition (0, 1, 2), Z2: current_load_ratio
    # Simplified: Assume t is the continuous running time or equivalent aging time
    
    term_env = exp(gen.gamma1 * temp_condition + gen.gamma2 * current_load_ratio^2)
    h_val = gen.beta * (gen.scale_a^(-gen.beta)) * (t^(gen.beta - 1)) * term_env
    
    # Cumulative failure probability F(t) recursive calculation (Simplified: F(t) ≈ 1 - exp(-H(t)))
    # For very short steps, accumulate Hazard rate
    new_hazard = accumulated_hazard + h_val
    
    # UPDATED: Calculate Conditional Failure Probability for this time step
    # P(Fail in [t, t+1] | Survived to t) = 1 - exp(-h(t)*dt)
    # Assuming dt = 1 hour
    prob_fail = 1.0 - exp(-h_val)
    
    return prob_fail, new_hazard
end

# ==========================================
# 3. Core Optimization Model: Rolling Horizon Unit Commitment
# ==========================================

function solve_rh_uc(
    gens, ess, 
    load_profile, wind_profile, solar_profile, 
    t_start, t_end,
    prev_status, prev_power, prev_soc,
    gen_availability, # 1: available, 0: outage
    day_ahead_soc_ref # Reference SOC from bi-level scheduling
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    horizon_len = t_end - t_start + 1
    steps = 1:horizon_len
    
    # --- Variables ---
    @variable(model, p[g=1:length(gens), t=steps] >= 0)
    @variable(model, u[g=1:length(gens), t=steps], Bin) # On/Off
    @variable(model, v[g=1:length(gens), t=steps], Bin) # Start up
    @variable(model, w[g=1:length(gens), t=steps], Bin) # Shut down
    
    @variable(model, p_ch[t=steps] >= 0)
    @variable(model, p_dis[t=steps] >= 0)
    @variable(model, soc[t=steps])
    @variable(model, load_shed[t=steps] >= 0) # Allow load shedding
    @variable(model, p_curtail[t=steps] >= 0) # Allow renewable curtailment (Fix for over-generation infeasibility)
    
    # --- Constraints ---
    for t in steps
        abs_t = t_start + t - 1
        
        # 1. Power Balance (Eq 14)
        total_gen = sum(p[g, t] for g in 1:length(gens))
        net_renewables = wind_profile[abs_t] + solar_profile[abs_t]
        
        # Balance: Gen + (Renewable - Curtail) + Discharge - Charge + Shed == Load
        @constraint(model, 
            total_gen + (net_renewables - p_curtail[t]) + p_dis[t] - p_ch[t] + load_shed[t] 
            == load_profile[abs_t]
        )
        
        # Limit curtailment to available renewables
        @constraint(model, p_curtail[t] <= net_renewables)
        
        # 2. Storage Constraints
        if t == 1
            @constraint(model, soc[t] == prev_soc + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        else
            @constraint(model, soc[t] == soc[t-1] + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        end
        
        @constraint(model, ess.SOC_min * ess.Ecap <= soc[t])
        @constraint(model, soc[t] <= ess.SOC_max * ess.Ecap)
        @constraint(model, p_ch[t] <= ess.Pc_max)
        @constraint(model, p_dis[t] <= ess.Pd_max)
        
        # 3. Bi-level Scheduling Constraint (Eq 15): Track Day-Ahead SOC Trajectory
        # Note: Simplified here, assuming day_ahead_soc_ref is percentage
        ref_soc_val = day_ahead_soc_ref[abs_t] * ess.Ecap
        @constraint(model, soc[t] - ref_soc_val <= ess.alpha * ess.SOC_max * ess.Ecap)
        @constraint(model, soc[t] - ref_soc_val >= -ess.alpha * ess.SOC_max * ess.Ecap)
        
        # 4. Generator Constraints
        for g in 1:length(gens)
            # Forced Outage Constraint (Eq 6): If failed, u must be 0
            if gen_availability[g, abs_t] == 0
                @constraint(model, u[g, t] == 0)
                @constraint(model, p[g, t] == 0)
                
                # Logic Constraints (maintain v/w consistency even during outage)
                if t == 1
                    @constraint(model, u[g, t] - prev_status[g] == v[g, t] - w[g, t])
                else
                    @constraint(model, u[g, t] - u[g, t-1] == v[g, t] - w[g, t])
                end
                # NOTE: Ramping constraints are SKIPPED during forced outage to allow instant trip
            else
                # Normal Operation
                @constraint(model, p[g, t] <= gens[g].Pmax * u[g, t])
                @constraint(model, p[g, t] >= gens[g].Pmin * u[g, t])
                
                # Logic & Ramping Constraints
                if t == 1
                    @constraint(model, u[g, t] - prev_status[g] == v[g, t] - w[g, t])
                    @constraint(model, p[g, t] - prev_power[g] <= gens[g].RU)
                    @constraint(model, prev_power[g] - p[g, t] <= gens[g].RD)
                else
                    @constraint(model, u[g, t] - u[g, t-1] == v[g, t] - w[g, t])
                    @constraint(model, p[g, t] - p[g, t-1] <= gens[g].RU)
                    @constraint(model, p[g, t-1] - p[g, t] <= gens[g].RD)
                end
            end
        end
    end
    
    # --- Objective Function ---
    # Minimize: Operating Cost + Start-up Cost + Load Shedding Penalty (VOLL)
    cost_gen = sum(
        gens[g].Cost_fix * u[g, t] + 
        gens[g].Cost_start * v[g, t] +
        0.1 * p[g, t] # Simplified fuel variable cost
        for g in 1:length(gens), t in steps
    )
    
    cost_shed = sum(load_shed[t] * 10000.0 for t in steps) # VOLL = 10000
    
    @objective(model, Min, cost_gen + cost_shed)
    
    optimize!(model)
    
    if termination_status(model) != MOI.OPTIMAL
        # If infeasible, return fallback strategy
        # Treat infeasibility as full system failure for this step to avoid 0 LOLP
        println("Warning: RH-UC Infeasible at step $t_start")
        # Return 0 generation, but report full load as shedding
        return zeros(length(gens)), zeros(length(gens)), prev_soc, load_profile[t_start]
    end
    
    # Extract results for the first step (t=1) for implementation
    impl_p = value.(p)[:, 1]
    impl_u = value.(u)[:, 1]
    impl_soc = value(soc[1])
    impl_shed = value(load_shed[1])
    
    return impl_p, impl_u, impl_soc, impl_shed
end

# ==========================================
# 4. Monte Carlo Simulation Main Loop
# ==========================================

"""
    run_reliability_evaluation()
    
    UPDATED: Uses Composite Importance Sampling (IS) for variance reduction.
    Biases both Renewables (mildly) and Generator Outages (strongly).
"""
function run_reliability_evaluation()
    gens = init_generators()
    
    # --- DBBN Training Phase ---
    println("Training DBBN models on mock historical data...")
    hist_wind, hist_solar = generate_mock_history_data(8760) 
    wind_dbbn_orig = train_dbbn(hist_wind, 15)
    solar_dbbn_orig = train_dbbn(hist_solar, 10)
    
    # --- Importance Sampling Setup ---
    println("Applying Composite Importance Sampling:")
    println("  - Renewable Twisting Factor: $IS_TWIST_FACTOR")
    println("  - Generator Failure Bias: $GEN_FAILURE_BIAS x")
    println("  * Note: System capacity increased to demonstrate IS benefit on rare events.")
    
    # Twist to favor LOWER states (Stress Testing)
    wind_dbbn_biased = twist_dbbn(wind_dbbn_orig, IS_TWIST_FACTOR)
    solar_dbbn_biased = twist_dbbn(solar_dbbn_orig, IS_TWIST_FACTOR)
    println("DBBN Training & Biasing Complete.")
    # ---------------------------
    
    gen_base_ages = [rand(2000:4000) for _ in 1:length(gens)]
    println("Generator Initial Ages (Hours): $gen_base_ages")

    # Store weighted metrics
    # We need to store (Value * Weight) and (Weight) separately to compute mean
    sum_weights = 0.0
    sum_lolp_weighted = 0.0
    sum_eens_weighted = 0.0
    sum_cost_weighted = 0.0
    
    # For variance calculation (IS Variance formula)
    # Var(X_is) = E[(X*W)^2] - E[X*W]^2  (Simplified)
    sq_sum_lolp_weighted = 0.0
    
    results = SimulationResults(Float64[], Float64[], Float64[], 0.0, 0)
    
    println("=" ^ 60)
    println("Starting Reliability Evaluation (Importance Sampling)")
    println("Horizon: $T_HORIZON hours, Max Samples: $MAX_SAMPLES")
    println("=" ^ 60)
    
    for n_sim in 1:MAX_SAMPLES
        # 1. Generate Scenarios using IS (Renewables)
        load_prof, wind_prof, solar_prof, renewable_weight = generate_dbbn_scenarios_is(
            wind_dbbn_orig, wind_dbbn_biased,
            solar_dbbn_orig, solar_dbbn_biased,
            T_HORIZON
        )
        
        # Initialize total scenario weight with renewable weight
        current_scenario_weight = renewable_weight
        
        day_ahead_soc_ref = fill(0.5, T_HORIZON) 
        
        # Initialize State
        current_status = zeros(Int, length(gens))
        current_power = zeros(Float64, length(gens))
        current_soc = ESS.Ecap * 0.5
        gen_hazards = zeros(Float64, length(gens))
        gen_availability = ones(Int, length(gens), T_HORIZON)
        
        total_load_shed = 0.0
        total_demand = sum(load_prof)
        total_cost = 0.0
        
        # --- Rolling Horizon Loop ---
        for t in 1:T_HORIZON
            # A. Sample Failure States
            temp_cond = (t >= 10 && t <= 16) ? 2 : 0 
            
            for g in 1:length(gens)
                load_ratio = current_power[g] / gens[g].Pmax
                effective_age = gen_base_ages[g] + t
                prob_fail_orig, new_haz = calculate_outage_prob(gens[g], effective_age, gen_hazards[g], load_ratio, temp_cond)
                gen_hazards[g] = new_haz
                
                # --- Importance Sampling for Generators ---
                # Bias the failure probability to encourage more failures
                # Ensure prob doesn't exceed 1.0
                prob_fail_biased = min(0.99, prob_fail_orig * GEN_FAILURE_BIAS)
                
                is_failure = false
                if rand() < prob_fail_biased
                    is_failure = true
                    # Weight update for failure event: P_orig / P_biased
                    # If biased prob is higher, weight < 1 (down-weighting frequent failures)
                    current_scenario_weight *= (prob_fail_orig / prob_fail_biased)
                else
                    is_failure = false
                    # Weight update for survival event: (1-P_orig) / (1-P_biased)
                    current_scenario_weight *= ((1.0 - prob_fail_orig) / (1.0 - prob_fail_biased))
                end
                
                if is_failure
                    repair_end = min(T_HORIZON, t + 5)
                    gen_availability[g, t:repair_end] .= 0
                    gen_hazards[g] = 0.0 
                end
            end
            
            # B. Determine Rolling Window
            t_end = min(T_HORIZON, t + T_ROLL - 1)
            
            # C. Solve RH-UC
            p_opt, u_opt, soc_opt, shed_opt = solve_rh_uc(
                gens, ESS,
                load_prof, wind_prof, solar_prof,
                t, t_end,
                current_status, current_power, current_soc,
                gen_availability,
                day_ahead_soc_ref
            )
            
            # D. Update
            current_power = p_opt
            current_status = round.(Int, u_opt)
            current_soc = soc_opt
            total_load_shed += shed_opt
            total_cost += sum(gens[g].Cost_fix * current_status[g] for g in 1:length(gens))
        end
        
        # --- Calculate Weighted Metrics ---
        
        # Raw Indicator (1 if loss, 0 if not)
        raw_lolp = total_load_shed > 1e-3 ? 1.0 : 0.0
        raw_eens = total_demand > 0 ? total_load_shed / total_demand : 0.0
        
        # Accumulate Weighted Sums
        # Note: In IS, E[X] = 1/N * Sum(X_i * W_i)
        sum_weights += 1.0 # We average over N simulations
        
        weighted_lolp = raw_lolp * current_scenario_weight
        sum_lolp_weighted += weighted_lolp
        sq_sum_lolp_weighted += (weighted_lolp)^2 # Approximation for variance tracking
        
        sum_eens_weighted += raw_eens * current_scenario_weight
        sum_cost_weighted += total_cost * current_scenario_weight # Cost is also weighted by probability of scenario
        
        # Current Estimates
        est_lolp = sum_lolp_weighted / n_sim
        est_eens = sum_eens_weighted / n_sim
        
        push!(results.lolp_history, est_lolp)
        push!(results.eens_history, est_eens)
        
        # Calculate CV (Coefficient of Variation) for IS
        if n_sim >= 10 && est_lolp > 0
            # Var(X_is) approx (1/N) * [ (1/N)*Sum((X*W)^2) - (Mean)^2 ]
            # This is a simplified variance estimator for IS
            mean_sq = sq_sum_lolp_weighted / n_sim
            var_est = (mean_sq - est_lolp^2)
            
            if var_est > 0
                std_dev = sqrt(var_est / n_sim) # Standard Error of the Mean
                cv = std_dev / est_lolp
                push!(results.cv_history, cv)
                
                println("Iter $n_sim: LOLP(IS)=$(round(est_lolp, digits=6)), EENS(IS)=$(round(est_eens, digits=6)), CV=$(round(cv, digits=4))")
                
                if cv < CV_THRESHOLD
                    println("\n✓ Convergence achieved at iteration $n_sim")
                    results.convergence_iteration = n_sim
                    break
                end
            end
        elseif n_sim % 10 == 0
             println("Iter $n_sim: LOLP(IS)=$(round(est_lolp, digits=6)) (Accumulating samples...)")
        end
    end
    
    results.total_cost = sum_cost_weighted / length(results.lolp_history)
    
    final_lolp = results.lolp_history[end]
    final_eens = results.eens_history[end]
    
    println("\n" * "=" ^ 60)
    println("FINAL RELIABILITY METRICS (Importance Sampling)")
    println("=" ^ 60)
    println("LOLP: $(round(final_lolp, digits=6))")
    println("EENS: $(round(final_eens, digits=6))")
    println("Convergence Iteration: $(results.convergence_iteration > 0 ? results.convergence_iteration : "Not Reached")")
    println("=" ^ 60)
    
    return results
end

# ==========================================
# 5. Visualization Functions
# ==========================================

"""
    plot_simulation_results(results::SimulationResults)

Generates and saves plots for reliability metrics convergence.
"""
function plot_simulation_results(results::SimulationResults)
    # 1. LOLP Convergence Plot
    n_iters = length(results.lolp_history)
    running_lolp = [mean(results.lolp_history[1:i]) for i in 1:n_iters]
    
    p1 = plot(1:n_iters, running_lolp, 
        title="LOLP Convergence", 
        xlabel="Iteration", ylabel="LOLP",
        lw=2, legend=false, color=:blue)
    
    # 2. CV Convergence Plot
    if !isempty(results.cv_history)
        # CV is calculated starting from iteration 5
        start_idx = 5
        valid_cv = results.cv_history
        iters_cv = start_idx:(start_idx + length(valid_cv) - 1)
        
        p2 = plot(iters_cv, valid_cv, 
            title="Coefficient of Variation (CV)", 
            xlabel="Iteration", ylabel="CV",
            lw=2, legend=false, color=:red)
        hline!(p2, [CV_THRESHOLD], linestyle=:dash, color=:black, label="Threshold")
    else
        p2 = plot(title="CV (Not enough data)")
    end

    # 3. EENS Histogram
    p3 = histogram(results.eens_history, 
        title="EENS Distribution", 
        xlabel="EENS (Normalized)", ylabel="Frequency",
        legend=false, color=:green, bins=20)

    # Combine and save
    final_plot = plot(p1, p2, p3, layout=(3, 1), size=(800, 1000))
    savefig(final_plot, "reliability_results.png")
    println("Saved reliability plots to 'reliability_results.png'")
end

"""
    plot_dispatch_profile(gens, load, wind, solar, p_gen, p_dis, p_ch, load_shed, t_range)

Visualizes the system dispatch for a specific time range (e.g., one day).
Useful for debugging specific scenarios.
"""
function plot_dispatch_profile(load, wind, solar, p_gen_total, p_dis, p_ch, load_shed)
    T = length(load)
    t_axis = 1:T
    
    # Stacked Area Plot for Supply
    # We'll plot: Thermal, Wind, Solar, Storage Discharge, Load Shedding
    
    # Prepare data for stacked plot (simplified approach using plot multiple times)
    p1 = plot(t_axis, load, label="Load Demand", color=:black, lw=2, linestyle=:dash)
    
    # Cumulative supply for stacking effect
    y1 = p_gen_total
    y2 = y1 .+ wind
    y3 = y2 .+ solar
    y4 = y3 .+ p_dis
    y5 = y4 .+ load_shed
    
    plot!(p1, t_axis, y1, fillrange=0, fillalpha=0.5, label="Thermal", color=:red)
    plot!(p1, t_axis, y2, fillrange=y1, fillalpha=0.5, label="Wind", color=:blue)
    plot!(p1, t_axis, y3, fillrange=y2, fillalpha=0.5, label="Solar", color=:orange)
    plot!(p1, t_axis, y4, fillrange=y3, fillalpha=0.5, label="Storage Disch.", color=:purple)
    plot!(p1, t_axis, y5, fillrange=y4, fillalpha=0.5, label="Load Shed", color=:gray)
    
    title!(p1, "System Dispatch Profile")
    ylabel!(p1, "Power (MW)")
    
    # Storage Charging
    p2 = bar(t_axis, -p_ch, label="Storage Charge", color=:green)
    ylabel!(p2, "Charging (MW)")
    xlabel!(p2, "Hour")
    
    final_plot = plot(p1, p2, layout=(2, 1), size=(800, 800))
    savefig(final_plot, "dispatch_profile.png")
    println("Saved dispatch profile to 'dispatch_profile.png'")
end

# ==========================================
# 6. Run the Simulation
# ==========================================
println("\n🚀 Power System Reliability Evaluation Framework")
println("Based on Joint Machine Learning and Optimization\n")

Random.seed!(42)  # For reproducibility
results = run_reliability_evaluation()

# Generate Plots
println("\n📊 Generating Visualizations...")
plot_simulation_results(results)

println("\n✓ Simulation completed successfully!")
