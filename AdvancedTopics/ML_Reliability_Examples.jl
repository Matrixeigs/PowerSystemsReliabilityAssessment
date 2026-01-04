"""
Machine Learning Enhanced Reliability Evaluation - Complete Implementation
===========================================================================

This file contains comprehensive Julia implementations combining educational
clarity with production-ready features from the lecture and main.jl:

BASIC FEATURES (Educational):
1. DBBN-based renewable sampling
2. Proportional Hazard Model (PHM) for aging-aware failures
3. Rolling-Horizon Unit Commitment
4. Sequential Monte Carlo Integration

ADVANCED FEATURES (Production):
5. Full MILP UC formulation (with binary variables)
6. Composite Importance Sampling for variance reduction
7. DBBN training from historical data
8. Emergency slack variables for feasibility
9. Bi-level storage dispatch

Author: Based on Lin & Shang (IEEE TPWRS 2021)
Course: Electric Power System Reliability Assessment
"""

using Random, Distributions, Statistics, JuMP, HiGHS, Plots

# Configuration: Set to true to enable advanced features
const USE_MILP = true              # Use full MILP formulation (vs linear relaxation)
const USE_IMPORTANCE_SAMPLING = false  # Use IS for variance reduction
const USE_DBBN_TRAINING = true     # Train DBBN from mock data (vs simple sampling)

println("="^70)
println("ML-Enhanced Reliability Evaluation - Configuration")
println("="^70)
println("MILP Formulation:      $(USE_MILP ? "✓ Enabled" : "✗ Disabled (Linear Relaxation)")")
println("Importance Sampling:   $(USE_IMPORTANCE_SAMPLING ? "✓ Enabled" : "✗ Disabled (Standard MC)")")
println("DBBN Training:         $(USE_DBBN_TRAINING ? "✓ Enabled" : "✗ Disabled (Simple Model)")")
println("="^70)
println()

# =============================================================================
# PART 1: DBBN - Dynamic Bayesian Belief Network for Renewables
# =============================================================================

"""
    DiscreteRenewableDistribution

Represents a discrete probability distribution for renewable output at time t.
This is the OUTPUT of a trained DBBN model (from lecture slides Section 1).
"""
struct DiscreteRenewableDistribution
    bins::Vector{Float64}        # Discretized output levels (MW)
    probabilities::Vector{Float64}  # P(output = bins[i])
end

"""
    DiscreteDBBN

Advanced DBBN implementation with temporal dependencies.
Stores Conditional Probability Tables (CPT) as transition matrix.
"""
struct DiscreteDBBN
    state_centers::Vector{Float64}      # Physical value for each state
    transition_matrix::Matrix{Float64}  # P(State_t | State_{t-1})
    initial_dist::Vector{Float64}       # Initial state probability
end

# --- Basic DBBN Functions ---

"""
    sample_renewable(dist::DiscreteRenewableDistribution)

Inverse transform sampling from discrete distribution (Lecture Section 3.2).
"""
function sample_renewable(dist::DiscreteRenewableDistribution)
    # Build CDF
    cdf = cumsum(dist.probabilities)

    # Sample u ~ Uniform(0,1)
    u = rand()

    # Invert CDF: find bin where u falls
    idx = findfirst(c -> c >= u, cdf)

    return dist.bins[idx]
end

"""
    create_example_wind_distribution(hour::Int)

Creates a simple time-varying wind distribution.
In practice, this would come from a trained DBBN model.
"""
function create_example_wind_distribution(hour::Int)
    # Time-dependent mean (higher wind at night/early morning)
    time_factor = 0.7 + 0.3 * cos(2π * (hour - 3) / 24)

    # Discrete bins from 0 to 100 MW
    bins = collect(0:10:100)
    n_bins = length(bins)

    # Beta distribution parameters (shape changes with time)
    α = 2.0 + time_factor * 3.0
    β = 3.0 - time_factor * 1.5

    # Generate probabilities
    probs = [pdf(Beta(α, β), (i-1)/(n_bins-1)) for i in 1:n_bins]
    probs ./= sum(probs)  # Normalize

    return DiscreteRenewableDistribution(bins, probs)
end

# --- Advanced DBBN Functions (from main.jl) ---

"""
    generate_mock_history_data(n_hours::Int)

Simulates historical wind/solar data for DBBN training.
Uses AR(1) process for wind, daily cycle for solar.
"""
function generate_mock_history_data(n_hours::Int)
    # Wind Speed (Weibull + Autocorrelation)
    wind_speed = zeros(n_hours)
    wind_speed[1] = 5.0
    for t in 2:n_hours
        noise = rand(Normal(0, 1.5))
        wind_speed[t] = 0.9 * wind_speed[t-1] + 0.1 * 7.0 + noise
        wind_speed[t] = clamp(wind_speed[t], 0.0, 25.0)
    end

    # Solar Irradiance (Daily Cycle + Cloud noise)
    solar_irrad = zeros(n_hours)
    for t in 1:n_hours
        hour = (t - 1) % 24
        if 6 <= hour <= 18
            max_solar = 1000.0 * sin((hour - 6) * π / 12)
            cloud_factor = rand(Beta(2, 1))
            solar_irrad[t] = max_solar * cloud_factor
        else
            solar_irrad[t] = 0.0
        end
    end

    return wind_speed, solar_irrad
end

"""
    discretize_data(data::Vector{Float64}, n_states::Int)

Discretizes continuous time series into bins for CPT learning.
"""
function discretize_data(data::Vector{Float64}, n_states::Int)
    min_val, max_val = minimum(data), maximum(data)
    edges = range(min_val, max_val, length=n_states+1)
    state_centers = [(edges[i] + edges[i+1])/2 for i in 1:n_states]

    discrete_seq = Int[]
    for x in data
        idx = searchsortedlast(edges, x)
        idx = clamp(idx, 1, n_states)
        push!(discrete_seq, idx)
    end

    return discrete_seq, state_centers
end

"""
    train_dbbn(data::Vector{Float64}, n_states::Int=10)

Trains a DBBN model from historical data using maximum likelihood estimation.
"""
function train_dbbn(data::Vector{Float64}, n_states::Int=10)
    seq, centers = discretize_data(data, n_states)
    counts = zeros(Int, n_states, n_states)
    state_counts = zeros(Int, n_states)

    # Count transitions
    for t in eachindex(seq)[2:end]
        prev_s = seq[t-1]
        curr_s = seq[t]
        counts[prev_s, curr_s] += 1
        state_counts[prev_s] += 1
    end
    state_counts[seq[end]] += 1

    # Estimate transition probabilities
    trans_mat = zeros(Float64, n_states, n_states)
    for i in 1:n_states
        row_sum = sum(counts[i, :])
        if row_sum > 0
            trans_mat[i, :] = counts[i, :] ./ row_sum
        else
            trans_mat[i, i] = 1.0
        end
    end

    # Initial distribution
    init_dist = state_counts ./ sum(state_counts)
    if sum(init_dist) == 0
        init_dist = fill(1.0/n_states, n_states)
    end

    return DiscreteDBBN(centers, trans_mat, init_dist)
end

"""
    sample_dbbn(model::DiscreteDBBN, n_steps::Int)

Samples a trajectory from trained DBBN model.
"""
function sample_dbbn(model::DiscreteDBBN, n_steps::Int)
    trajectory = zeros(Float64, n_steps)
    current_state = rand(Categorical(model.initial_dist))
    trajectory[1] = model.state_centers[current_state]

    for t in 2:n_steps
        probs = model.transition_matrix[current_state, :]
        probs = probs ./ sum(probs)
        next_state = rand(Categorical(probs))
        trajectory[t] = model.state_centers[next_state]
        current_state = next_state
    end

    return trajectory
end

"""
    wind_power_curve(speed_ms)

Converts wind speed to power using typical turbine curve.
"""
function wind_power_curve(speed_ms)
    v_in = 3.0
    v_rated = 12.0
    v_out = 25.0
    p_rated = 450.0  # MW

    if speed_ms < v_in || speed_ms > v_out
        return 0.0
    elseif speed_ms >= v_rated
        return p_rated
    else
        return p_rated * ((speed_ms - v_in) / (v_rated - v_in))^3
    end
end

# --- Importance Sampling (Advanced) ---

"""
    twist_dbbn(model::DiscreteDBBN, twisting_factor::Float64)

Applies exponential twisting for Importance Sampling.
Positive factor favors lower states (stress testing).
"""
function twist_dbbn(model::DiscreteDBBN, twisting_factor::Float64)
    n_states = length(model.state_centers)
    new_trans = zeros(Float64, n_states, n_states)

    for i in 1:n_states
        denom = 0.0
        for j in 1:n_states
            bias = exp(-twisting_factor * (j / n_states))
            new_trans[i, j] = model.transition_matrix[i, j] * bias
            denom += new_trans[i, j]
        end
        if denom > 0
            new_trans[i, :] ./= denom
        else
            new_trans[i, i] = 1.0
        end
    end

    # Twist initial distribution
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
    sample_dbbn_with_likelihood(original, biased, n_steps)

Samples from biased DBBN and returns trajectory + likelihood ratio weight.
"""
function sample_dbbn_with_likelihood(original::DiscreteDBBN, biased::DiscreteDBBN, n_steps::Int)
    trajectory = zeros(Float64, n_steps)
    log_likelihood_ratio = 0.0

    # Sample initial state from biased
    current_state = rand(Categorical(biased.initial_dist))
    trajectory[1] = biased.state_centers[current_state]

    # Weight for initial state
    p_orig = original.initial_dist[current_state]
    p_bias = biased.initial_dist[current_state]
    if p_bias > 0 && p_orig > 0
        log_likelihood_ratio += log(p_orig) - log(p_bias)
    end

    # Sample transitions
    for t in 2:n_steps
        probs_bias = biased.transition_matrix[current_state, :]
        next_state = rand(Categorical(probs_bias))

        p_trans_orig = original.transition_matrix[current_state, next_state]
        p_trans_bias = biased.transition_matrix[current_state, next_state]

        if p_trans_bias > 0 && p_trans_orig > 0
            log_likelihood_ratio += log(p_trans_orig) - log(p_trans_bias)
        end

        trajectory[t] = biased.state_centers[next_state]
        current_state = next_state
    end

    return trajectory, exp(log_likelihood_ratio)
end

# =============================================================================
# PART 2: Proportional Hazard Model (PHM)
# =============================================================================

"""
    PHMParameters

Parameters for the Proportional Hazard Model (Lecture Section 2.3).
h(t, Z) = a^(-b) * b * t^(b-1) * exp(γ₁*Z₁ + γ₂*Z₂²)
"""
struct PHMParameters
    a::Float64      # Weibull scale parameter (hours)
    b::Float64      # Weibull shape parameter (aging effect)
    γ₁::Float64     # Temperature sensitivity coefficient
    γ₂::Float64     # Loading sensitivity coefficient
end

"""
    compute_hazard_rate(phm, t, temp_category, loading_rate)

Computes instantaneous hazard rate using PHM (Lecture Eq. 2.3).
"""
function compute_hazard_rate(phm::PHMParameters, t::Float64, temp_category::Int, loading_rate::Float64)
    # Baseline Weibull hazard
    baseline = (phm.a^(-phm.b)) * phm.b * (t^(phm.b - 1))

    # Stress multiplier
    stress_multiplier = exp(phm.γ₁ * temp_category + phm.γ₂ * loading_rate^2)

    return baseline * stress_multiplier
end

"""
    update_outage_probability(F_prev, phm, t, temp, loading)

Recursive outage probability update (Lecture Eq. 2.4).
F(t) = 1 - exp(-∫h dx) + F(t-1)*exp(-∫h dx)
"""
function update_outage_probability(F_prev::Float64, phm::PHMParameters,
                                  t::Float64, temp_category::Int, loading_rate::Float64)
    h = compute_hazard_rate(phm, t, temp_category, loading_rate)

    # Δt = 1 hour assumption
    Δt = 1.0
    survival_prob = exp(-h * Δt)

    # Recursive formula
    F_new = 1.0 - survival_prob * (1.0 - F_prev)

    return F_new, h
end

# =============================================================================
# PART 3: Generator and System Models
# =============================================================================

"""
    Generator

Generator struct with complete UC and aging parameters.
Supports both linear relaxation and full MILP formulation.
"""
struct Generator
    id::Int
    Pmax::Float64
    Pmin::Float64
    RD::Float64        # Ramp Down (MW/h)
    RU::Float64        # Ramp Up (MW/h)
    TU::Int            # Min Up Time (hours) - only for MILP
    TD::Int            # Min Down Time (hours) - only for MILP
    Cost_fix::Float64  # Fixed cost ($/h)
    Cost_start::Float64  # Startup cost ($)
    # PHM parameters
    beta::Float64
    scale_a::Float64
    gamma1::Float64
    gamma2::Float64
end

"""
    Storage

Energy Storage System parameters (Lecture Section 5).
"""
struct Storage
    Ecap::Float64       # Energy capacity (MWh)
    Pc_max::Float64     # Charging power limit (MW)
    Pd_max::Float64     # Discharging power limit (MW)
    eta_c::Float64      # Charging efficiency
    eta_d::Float64      # Discharging efficiency
    SOC_min::Float64    # Minimum SOC (fraction)
    SOC_max::Float64    # Maximum SOC (fraction)
    alpha::Float64      # RH tracking tolerance
end

"""
    init_generators()

Initialize generator fleet with realistic parameters.
Total capacity ~1150 MW for positive reserve margin.
"""
function init_generators()
    return [
        Generator(1, 600.0, 150.0, 200.0, 200.0, 8, 8, 400.0, 5000.0, 3.0, 5000.0, 0.5, 1.0),
        Generator(2, 350.0, 50.0, 100.0, 100.0, 4, 4, 200.0, 2000.0, 3.0, 4000.0, 0.5, 1.0),
        Generator(3, 200.0, 20.0, 60.0, 60.0, 1, 1, 100.0, 500.0, 3.0, 3000.0, 0.5, 1.0)
    ]
end

# System Parameters
const ESS = Storage(300.0, 60.0, 60.0, 0.95, 0.95, 0.1, 0.9, 1.0)
const T_HORIZON = 24
const T_ROLL = 6
const LOAD_PEAK = 1000.0
const CV_THRESHOLD = 0.05
const MAX_SAMPLES = 2000

# =============================================================================
# PART 4: Rolling-Horizon Unit Commitment
# =============================================================================

"""
    solve_rh_uc_linear(...)

Simplified linear relaxation UC (for educational purposes).
No binary variables - generators can operate at any level 0 to Pmax.
"""
function solve_rh_uc_linear(
    gens, ess,
    load_profile, wind_profile, solar_profile,
    t_start, t_end,
    prev_power, prev_soc,
    gen_availability,
    day_ahead_soc_ref
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    horizon_len = t_end - t_start + 1
    steps = 1:horizon_len
    n_gen = length(gens)

    # Variables
    @variable(model, p[1:n_gen, t=steps] >= 0)
    @variable(model, p_ch[t=steps] >= 0)
    @variable(model, p_dis[t=steps] >= 0)
    @variable(model, soc[t=steps])
    @variable(model, load_shed[t=steps] >= 0)
    @variable(model, p_curtail[t=steps] >= 0)

    # Constraints
    for t in steps
        abs_t = t_start + t - 1

        # Power Balance
        total_gen = sum(p[g, t] for g in 1:n_gen)
        net_renewables = wind_profile[abs_t] + solar_profile[abs_t]
        @constraint(model,
            total_gen + (net_renewables - p_curtail[t]) + p_dis[t] - p_ch[t] + load_shed[t]
            == load_profile[abs_t]
        )
        @constraint(model, p_curtail[t] <= net_renewables)

        # Storage
        if t == 1
            @constraint(model, soc[t] == prev_soc + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        else
            @constraint(model, soc[t] == soc[t-1] + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        end
        @constraint(model, ess.SOC_min * ess.Ecap <= soc[t] <= ess.SOC_max * ess.Ecap)
        @constraint(model, p_ch[t] <= ess.Pc_max)
        @constraint(model, p_dis[t] <= ess.Pd_max)

        # SOC Tracking
        ref_soc_val = day_ahead_soc_ref[abs_t] * ess.Ecap
        @constraint(model, soc[t] - ref_soc_val <= ess.alpha * ess.SOC_max * ess.Ecap)
        @constraint(model, soc[t] - ref_soc_val >= -ess.alpha * ess.SOC_max * ess.Ecap)

        # Generator Limits
        for g in 1:n_gen
            if gen_availability[g, abs_t] == 0
                @constraint(model, p[g, t] == 0)
            else
                @constraint(model, p[g, t] <= gens[g].Pmax)
                # Note: Pmin removed for linear relaxation feasibility
            end
        end

        # Ramp Constraints
        for g in 1:n_gen
            if gen_availability[g, abs_t] == 1
                if t == 1 && prev_power[g] > 0
                    @constraint(model, p[g, t] - prev_power[g] <= gens[g].RU)
                    @constraint(model, prev_power[g] - p[g, t] <= gens[g].RD)
                elseif t > 1
                    @constraint(model, p[g, t] - p[g, t-1] <= gens[g].RU)
                    @constraint(model, p[g, t-1] - p[g, t] <= gens[g].RD)
                end
            end
        end
    end

    # Objective
    cost_gen = sum(0.1 * p[g, t] for g in 1:n_gen, t in steps)
    cost_shed = sum(load_shed[t] * 10000.0 for t in steps)
    @objective(model, Min, cost_gen + cost_shed)

    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        return zeros(n_gen), prev_soc, load_profile[t_start]
    end

    return value.(p)[:, 1], value(soc[1]), value(load_shed[1])
end

"""
    solve_rh_uc_milp(...)

Full MILP formulation with binary variables and emergency slack.
Production-ready implementation from main.jl.
"""
function solve_rh_uc_milp(
    gens, ess,
    load_profile, wind_profile, solar_profile,
    t_start, t_end,
    prev_status, prev_power, prev_soc,
    gen_availability,
    day_ahead_soc_ref
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    horizon_len = t_end - t_start + 1
    steps = 1:horizon_len
    n_gen = length(gens)

    # Variables
    @variable(model, p[1:n_gen, t=steps] >= 0)
    @variable(model, u[1:n_gen, t=steps], Bin)  # On/Off
    @variable(model, v[1:n_gen, t=steps], Bin)  # Start up
    @variable(model, w[1:n_gen, t=steps], Bin)  # Shut down
    @variable(model, p_ch[t=steps] >= 0)
    @variable(model, p_dis[t=steps] >= 0)
    @variable(model, soc[t=steps])
    @variable(model, load_shed[t=steps] >= 0)
    @variable(model, p_curtail[t=steps] >= 0)
    @variable(model, pmin_slack[1:n_gen, t=steps] >= 0)  # Emergency slack

    # Constraints
    for t in steps
        abs_t = t_start + t - 1

        # Power Balance
        total_gen = sum(p[g, t] for g in 1:n_gen)
        net_renewables = wind_profile[abs_t] + solar_profile[abs_t]
        @constraint(model,
            total_gen + (net_renewables - p_curtail[t]) + p_dis[t] - p_ch[t] + load_shed[t]
            == load_profile[abs_t]
        )
        @constraint(model, p_curtail[t] <= net_renewables)

        # Storage (same as linear)
        if t == 1
            @constraint(model, soc[t] == prev_soc + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        else
            @constraint(model, soc[t] == soc[t-1] + p_ch[t]*ess.eta_c - p_dis[t]/ess.eta_d)
        end
        @constraint(model, ess.SOC_min * ess.Ecap <= soc[t] <= ess.SOC_max * ess.Ecap)
        @constraint(model, p_ch[t] <= ess.Pc_max)
        @constraint(model, p_dis[t] <= ess.Pd_max)

        # SOC Tracking
        ref_soc_val = day_ahead_soc_ref[abs_t] * ess.Ecap
        @constraint(model, soc[t] - ref_soc_val <= ess.alpha * ess.SOC_max * ess.Ecap)
        @constraint(model, soc[t] - ref_soc_val >= -ess.alpha * ess.SOC_max * ess.Ecap)

        # Generator Constraints
        for g in 1:n_gen
            if gen_availability[g, abs_t] == 0
                # Forced Outage
                @constraint(model, u[g, t] == 0)
                @constraint(model, p[g, t] == 0)

                # Logic (maintain consistency)
                if t == 1
                    @constraint(model, u[g, t] - prev_status[g] == v[g, t] - w[g, t])
                else
                    @constraint(model, u[g, t] - u[g, t-1] == v[g, t] - w[g, t])
                end
            else
                # Normal Operation
                @constraint(model, p[g, t] <= gens[g].Pmax * u[g, t])
                # Pmin with emergency slack (heavily penalized)
                @constraint(model, p[g, t] + pmin_slack[g, t] >= gens[g].Pmin * u[g, t])

                # Logic & Ramping
                if t == 1
                    @constraint(model, u[g, t] - prev_status[g] == v[g, t] - w[g, t])
                    if prev_status[g] > 0
                        # Relax ramp during transitions
                        @constraint(model, p[g, t] - prev_power[g] <= gens[g].RU + gens[g].Pmax * v[g, t])
                        @constraint(model, prev_power[g] - p[g, t] <= gens[g].RD + gens[g].Pmax * w[g, t])
                    end
                else
                    @constraint(model, u[g, t] - u[g, t-1] == v[g, t] - w[g, t])
                    @constraint(model, p[g, t] - p[g, t-1] <= gens[g].RU + gens[g].Pmax * v[g, t])
                    @constraint(model, p[g, t-1] - p[g, t] <= gens[g].RD + gens[g].Pmax * w[g, t])
                end
            end
        end
    end

    # Objective with slack penalty
    cost_gen = sum(
        gens[g].Cost_fix * u[g, t] +
        gens[g].Cost_start * v[g, t] +
        0.1 * p[g, t]
        for g in 1:n_gen, t in steps
    )
    cost_shed = sum(load_shed[t] * 10000.0 for t in steps)
    cost_slack = sum(pmin_slack[g, t] * 5000.0 for g in 1:n_gen, t in steps)
    @objective(model, Min, cost_gen + cost_shed + cost_slack)

    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        status = termination_status(model)
        println("⚠️  Warning: MILP UC $status at step $t_start")
        println("   Load: $(round(load_profile[t_start], digits=1)) MW")
        println("   Renewables: Wind=$(round(wind_profile[t_start], digits=1)), Solar=$(round(solar_profile[t_start], digits=1)) MW")
        println("   Available Gens: $(sum(gen_availability[:, t_start]))/$(n_gen)")
        return zeros(n_gen), zeros(Int, n_gen), prev_soc, load_profile[t_start]
    end

    return value.(p)[:, 1], round.(Int, value.(u)[:, 1]), value(soc[1]), value(load_shed[1])
end

# Wrapper function to select UC formulation
function solve_rh_uc(args...)
    if USE_MILP
        return solve_rh_uc_milp(args...)
    else
        p, soc, shed = solve_rh_uc_linear(args...)
        u = zeros(Int, length(args[1]))  # Dummy status for linear relaxation
        return p, u, soc, shed
    end
end

# =============================================================================
# PART 5: Sequential Monte Carlo Simulation
# =============================================================================

"""
    SimulationResults

Container for simulation outputs and convergence tracking.
"""
mutable struct SimulationResults
    lolp_history::Vector{Float64}
    cv_history::Vector{Float64}
    eens_history::Vector{Float64}
    total_cost::Float64
    convergence_iteration::Int
end

"""
    run_reliability_evaluation()

Main Sequential Monte Carlo loop integrating all components.
Supports both standard MC and Importance Sampling.
"""
function run_reliability_evaluation()
    gens = init_generators()
    gen_base_ages = [rand(2000:4000) for _ in 1:length(gens)]

    # Initialize DBBN models
    if USE_DBBN_TRAINING
        println("Training DBBN models from mock historical data...")
        hist_wind, hist_solar = generate_mock_history_data(8760)
        wind_dbbn = train_dbbn(hist_wind, 15)
        solar_dbbn = train_dbbn(hist_solar, 10)

        if USE_IMPORTANCE_SAMPLING
            println("Applying Importance Sampling (twisting factor: 0.05)")
            wind_dbbn_biased = twist_dbbn(wind_dbbn, 0.05)
            solar_dbbn_biased = twist_dbbn(solar_dbbn, 0.05)
        end
        println("DBBN setup complete.\n")
    end

    # Simulation tracking
    sum_lolp_weighted = 0.0
    sum_eens_weighted = 0.0
    sq_sum_lolp_weighted = 0.0

    results = SimulationResults(Float64[], Float64[], Float64[], 0.0, 0)

    println("="^60)
    println("Starting Reliability Evaluation")
    println("Horizon: $T_HORIZON hours, Max Samples: $MAX_SAMPLES")
    println("="^60)

    for n_sim in 1:MAX_SAMPLES
        # Generate scenario
        if USE_DBBN_TRAINING && USE_IMPORTANCE_SAMPLING
            wind_speed, w_wind = sample_dbbn_with_likelihood(wind_dbbn, wind_dbbn_biased, T_HORIZON)
            solar_irrad, w_solar = sample_dbbn_with_likelihood(solar_dbbn, solar_dbbn_biased, T_HORIZON)
            scenario_weight = w_wind * w_solar
        elseif USE_DBBN_TRAINING
            wind_speed = sample_dbbn(wind_dbbn, T_HORIZON)
            solar_irrad = sample_dbbn(solar_dbbn, T_HORIZON)
            scenario_weight = 1.0
        else
            # Simple sampling
            wind_speed = [sample_renewable(create_example_wind_distribution(h)) for h in 1:T_HORIZON]
            solar_irrad = zeros(T_HORIZON)
            scenario_weight = 1.0
        end

        # Convert to power
        wind_prof = USE_DBBN_TRAINING ? wind_power_curve.(wind_speed) : wind_speed
        solar_capacity = 200.0
        solar_prof = (solar_irrad ./ 1000.0) .* solar_capacity

        # Apply day/night cycle for solar
        for t in 1:T_HORIZON
            hour = (t - 1) % 24
            if hour < 6 || hour > 18
                solar_prof[t] = 0.0
            end
        end

        # Load profile
        load_prof = [LOAD_PEAK * (0.6 + 0.4 * sin((t-8)*π/12)) + rand(Normal(0, 20)) for t in 1:T_HORIZON]

        # Initialize system state
        current_status = zeros(Int, length(gens))
        current_power = zeros(Float64, length(gens))
        current_soc = ESS.Ecap * 0.5
        gen_F = zeros(Float64, length(gens))
        gen_availability = ones(Int, length(gens), T_HORIZON)

        total_load_shed = 0.0
        total_demand = sum(load_prof)

        # Day-ahead SOC reference
        day_ahead_soc_ref = fill(0.5, T_HORIZON)

        # Rolling Horizon Loop
        for t in 1:T_HORIZON
            # Temperature condition
            temp_cond = (t >= 10 && t <= 16) ? 2 : 0

            # Sample failures (with PHM)
            for g in eachindex(gens)
                if gen_availability[g, t] == 1
                    load_ratio = current_power[g] > 0 ? current_power[g] / gens[g].Pmax : 0.0
                    effective_age = Float64(gen_base_ages[g] + t)

                    # Update outage probability
                    gen_F[g], _ = update_outage_probability(
                        gen_F[g],
                        PHMParameters(gens[g].scale_a, gens[g].beta, gens[g].gamma1, gens[g].gamma2),
                        effective_age,
                        temp_cond,
                        load_ratio
                    )

                    # Sample failure
                    if rand() < gen_F[g]
                        repair_end = min(T_HORIZON, t + 5)
                        gen_availability[g, t:repair_end] .= 0
                        gen_F[g] = 0.0  # Reset after failure
                    end
                end
            end

            # Solve RH-UC
            t_end = min(T_HORIZON, t + T_ROLL - 1)

            if USE_MILP
                p_opt, u_opt, soc_opt, shed_opt = solve_rh_uc(
                    gens, ESS,
                    load_prof, wind_prof, solar_prof,
                    t, t_end,
                    current_status, current_power, current_soc,
                    gen_availability,
                    day_ahead_soc_ref
                )
                current_status = u_opt
            else
                p_opt, soc_opt, shed_opt = solve_rh_uc(
                    gens, ESS,
                    load_prof, wind_prof, solar_prof,
                    t, t_end,
                    current_power, current_soc,
                    gen_availability,
                    day_ahead_soc_ref
                )
            end

            current_power = p_opt
            current_soc = soc_opt
            total_load_shed += shed_opt
        end

        # Calculate metrics
        raw_lolp = total_load_shed > 1e-3 ? 1.0 : 0.0
        raw_eens = total_demand > 0 ? total_load_shed / total_demand : 0.0

        # Weighted accumulation
        weighted_lolp = raw_lolp * scenario_weight
        sum_lolp_weighted += weighted_lolp
        sum_eens_weighted += raw_eens * scenario_weight
        sq_sum_lolp_weighted += weighted_lolp^2

        # Current estimates
        est_lolp = sum_lolp_weighted / n_sim
        est_eens = sum_eens_weighted / n_sim

        push!(results.lolp_history, est_lolp)
        push!(results.eens_history, est_eens)

        # Convergence check
        if n_sim >= 10 && est_lolp > 0
            mean_sq = sq_sum_lolp_weighted / n_sim
            var_est = mean_sq - est_lolp^2

            if var_est > 0
                std_dev = sqrt(var_est / n_sim)
                cv = std_dev / est_lolp
                push!(results.cv_history, cv)

                if n_sim % 10 == 0
                    println("Iter $n_sim: LOLP=$(round(est_lolp, digits=6)), EENS=$(round(est_eens, digits=6)), CV=$(round(cv, digits=4))")
                end

                if cv < CV_THRESHOLD
                    println("\n✓ Convergence achieved at iteration $n_sim")
                    results.convergence_iteration = n_sim
                    break
                end
            end
        elseif n_sim % 20 == 0
            println("Iter $n_sim: LOLP=$(round(est_lolp, digits=6)), EENS=$(round(est_eens, digits=6))")
        end
    end

    println("\n" * "="^60)
    println("FINAL RELIABILITY METRICS")
    println("="^60)
    println("LOLP: $(round(results.lolp_history[end], digits=6))")
    println("EENS: $(round(results.eens_history[end], digits=6))")
    println("Convergence: $(results.convergence_iteration > 0 ? "Yes (iter $(results.convergence_iteration))" : "Not reached")")
    println("="^60)

    return results
end

"""
    plot_results(results::SimulationResults)

Generates convergence plots.
"""
function plot_results(results::SimulationResults)
    n_iters = length(results.lolp_history)

    p1 = plot(1:n_iters, results.lolp_history,
        title="LOLP Convergence",
        xlabel="Iteration", ylabel="LOLP",
        lw=2, legend=false, color=:blue)

    if !isempty(results.cv_history)
        start_idx = 10
        iters_cv = start_idx:(start_idx + length(results.cv_history) - 1)

        p2 = plot(iters_cv, results.cv_history,
            title="Coefficient of Variation",
            xlabel="Iteration", ylabel="CV",
            lw=2, legend=false, color=:red)
        hline!(p2, [CV_THRESHOLD], linestyle=:dash, color=:black, label="Target")
    else
        p2 = plot(title="CV (Not enough data)")
    end

    p3 = histogram(results.eens_history,
        title="EENS Distribution",
        xlabel="EENS", ylabel="Frequency",
        legend=false, color=:green, bins=20)

    final_plot = plot(p1, p2, p3, layout=(3, 1), size=(800, 1000))
    savefig(final_plot, "reliability_convergence.png")
    println("\n✓ Plot saved as 'reliability_convergence.png'")
end

# =============================================================================
# Main Execution
# =============================================================================

println("\n🚀 Power System Reliability Evaluation Framework")
println("Integrating DBBN + PHM + RH-UC + Sequential MCS\n")

Random.seed!(42)
results = run_reliability_evaluation()

println("\n📊 Generating Visualizations...")
plot_results(results)

println("\n✓ Simulation completed successfully!")
println("\nFor educational step-by-step examples, run: julia ML_Tutorial_StepByStep.jl")
println("For detailed guide, see: ML_Reliability_Guide.md")
