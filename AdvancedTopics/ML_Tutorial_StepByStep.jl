"""
Step-by-Step Tutorial: ML-Enhanced Reliability Evaluation
==========================================================

This tutorial walks through each concept individually with simple examples.
Perfect for first-time learners!

Each section can be run independently.
"""

using Random, Distributions, Plots

println("="^70)
println("ML-Enhanced Reliability Tutorial - Step by Step")
println("="^70)
println()

# =============================================================================
# STEP 1: Understanding DBBN - Sampling from Discrete Distributions
# =============================================================================

println("\n### STEP 1: DBBN - Discrete Renewable Sampling ###\n")

"""
Problem: Wind power is uncertain and time-dependent.
Traditional: Use Gaussian distribution (parametric)
DBBN: Use discrete probability distribution (non-parametric, captures tails better)
"""

# Simple example: Wind at 2 PM
println("Example: Wind power distribution at 2 PM (high wind likely)")
wind_bins = [0, 20, 40, 60, 80, 100]  # MW
wind_probs = [0.05, 0.1, 0.2, 0.35, 0.25, 0.05]  # Probabilities

println("Wind bins (MW):  ", wind_bins)
println("Probabilities:   ", wind_probs)
println("Sum of probs:    ", sum(wind_probs), " (must be 1.0)")

# Inverse Transform Sampling (Algorithm from slide 3.2)
function sample_wind_simple(bins, probs)
    # Step 1: Build CDF
    cdf = cumsum(probs)
    println("\nCDF: ", round.(cdf, digits=2))

    # Step 2: Sample u ~ Uniform(0,1)
    u = rand()
    println("Random u: ", round(u, digits=3))

    # Step 3: Find bin where u falls
    idx = findfirst(c -> c >= u, cdf)
    sampled_value = bins[idx]

    println("→ Sampled wind: $(sampled_value) MW (bin $idx)")
    return sampled_value
end

println("\n--- Single Sample ---")
sample_wind_simple(wind_bins, wind_probs)

println("\n--- Generate 1000 samples and check distribution ---")
samples = [sample_wind_simple(wind_bins, wind_probs) for _ in 1:1000]
println("Mean: ", round(mean(samples), digits=1), " MW")

# Visualize
histogram(samples, bins=wind_bins,
          xlabel="Wind Power (MW)",
          ylabel="Frequency",
          title="DBBN Samples (1000 samples)",
          legend=false)
savefig("step1_dbbn_samples.png")
println("\n✓ Plot saved as 'step1_dbbn_samples.png'")

println("\n" * "─"^70)
println("KEY INSIGHT: DBBN gives you the full distribution, not just mean!")
println("This matters for reliability - we care about tail events (low wind).")
println("─"^70)


# =============================================================================
# STEP 2: Understanding PHM - Stress Affects Failures
# =============================================================================

println("\n\n### STEP 2: PHM - Stress-Dependent Failure Rates ###\n")

"""
Problem: Generator failure rate should depend on operating stress
Traditional: λ = constant (e.g., 0.0001 failures/hour)
PHM: λ(t, stress) = f(age, temperature, loading)
"""

# Simple PHM example
function compute_failure_rate_traditional()
    return 0.0001  # failures/hour (constant)
end

function compute_failure_rate_phm(age_hours, temperature_category, loading_fraction)
    # Simplified PHM: h(t, Z) = baseline * stress_multiplier

    # Baseline (Weibull): increases with age
    baseline = 0.0001 * (age_hours / 8760)^0.2  # Slight aging effect

    # Stress multiplier
    temp_effect = exp(0.5 * temperature_category)  # 0=normal, 1=hot, 2=severe
    loading_effect = exp(2.0 * loading_fraction^2)  # Squared emphasizes overloading

    hazard = baseline * temp_effect * loading_effect
    return hazard
end

println("Compare failure rates:\n")

# Case 1: Young generator, light load, normal temperature
age1 = 1000.0  # hours
temp1 = 0      # normal
load1 = 0.5    # 50% loading
λ1 = compute_failure_rate_phm(age1, temp1, load1)
println("Case 1: Age=1000h, Temp=Normal, Load=50%")
println("  Traditional λ: $(compute_failure_rate_traditional())")
println("  PHM λ:        $(round(λ1, digits=6))")
println()

# Case 2: Old generator, heavy load, hot day
age2 = 50000.0  # hours
temp2 = 1       # hot
load2 = 0.95    # 95% loading (near max)
λ2 = compute_failure_rate_phm(age2, temp2, load2)
println("Case 2: Age=50000h, Temp=Hot, Load=95%")
println("  Traditional λ: $(compute_failure_rate_traditional())")
println("  PHM λ:        $(round(λ2, digits=6))")
println()

println("Ratio: Case 2 is $(round(λ2/λ1, digits=1))× more likely to fail!")

# Visualize loading effect
loadings = 0.1:0.05:1.0
hazards_normal = [compute_failure_rate_phm(10000.0, 0, l) for l in loadings]
hazards_hot = [compute_failure_rate_phm(10000.0, 1, l) for l in loadings]

plot(loadings, hazards_normal .* 1e6,
     label="Normal Temp",
     linewidth=2,
     xlabel="Loading Fraction",
     ylabel="Hazard Rate (×10⁻⁶)",
     title="PHM: Loading vs Failure Rate")
plot!(loadings, hazards_hot .* 1e6,
      label="Hot Temp",
      linewidth=2)
savefig("step2_phm_loading.png")
println("\n✓ Plot saved as 'step2_phm_loading.png'")

println("\n" * "─"^70)
println("KEY INSIGHT: Operational stress (loading, temperature) affects failures!")
println("This creates feedback: UC decisions → loading → failures → next UC")
println("─"^70)


# =============================================================================
# STEP 3: Understanding the Coupling Loop
# =============================================================================

println("\n\n### STEP 3: The Coupling Loop - UC ↔ PHM ↔ Reliability ###\n")

"""
Demonstrate the bidirectional feedback between UC dispatch and reliability.
"""

# Simplified scenario: 3 generators, varying load
struct SimpleGenerator
    id::Int
    max_power::Float64
    cost::Float64
end

generators = [
    SimpleGenerator(1, 100.0, 30.0),  # Cheap but will be loaded heavily
    SimpleGenerator(2, 80.0, 40.0),
    SimpleGenerator(3, 60.0, 50.0)
]

# Simple merit-order dispatch
function simple_dispatch(generators, load, failed_ids)
    # Sort by cost, exclude failed units
    available = filter(g -> !(g.id in failed_ids), generators)
    sort!(available, by=g -> g.cost)

    dispatch = Dict{Int, Float64}()
    remaining_load = load

    for gen in available
        power = min(gen.max_power, remaining_load)
        dispatch[gen.id] = power
        remaining_load -= power

        if remaining_load <= 0
            break
        end
    end

    load_shed = max(0.0, remaining_load)
    return dispatch, load_shed
end

println("Scenario: Load varies, observe loading and failures\n")

hours = 24
load_profile = [150 + 50*sin(2π*h/24) for h in 1:hours]  # Daily cycle

# Track state
generator_ages = [5000.0, 5000.0, 5000.0]  # hours since repair
failed_ids = Int[]
failure_events = []

for hour in 1:hours
    load = load_profile[hour]

    # Dispatch
    dispatch, ls = simple_dispatch(generators, load, failed_ids)

    # Calculate loading fractions
    loadings = Dict{Int, Float64}()
    for gen in generators
        if gen.id in failed_ids
            loadings[gen.id] = 0.0
        else
            loadings[gen.id] = get(dispatch, gen.id, 0.0) / gen.max_power
        end
    end

    # Update failure probabilities (simplified: check each hour)
    for gen in generators
        if gen.id in failed_ids
            continue  # Already failed
        end

        loading = loadings[gen.id]
        λ = compute_failure_rate_phm(generator_ages[gen.id], 0, loading)

        # Probability of failure in this 1-hour window: 1 - exp(-λ*1)
        prob_fail = 1 - exp(-λ * 1.0)

        if rand() < prob_fail
            push!(failed_ids, gen.id)
            push!(failure_events, (hour, gen.id, loading))
            println("  Hour $hour: Generator $(gen.id) FAILED (loading=$(round(loading*100, digits=0))%)")
        end
    end

    # Age all generators
    generator_ages .+= 1.0

    # Simple repair: failed units back after 4 hours
    # (In real code, track repair timers properly)
end

println("\nTotal failures: $(length(failure_events))")
if !isempty(failure_events)
    avg_loading_at_failure = mean([f[3] for f in failure_events])
    println("Average loading when failure occurred: $(round(avg_loading_at_failure*100, digits=0))%")
    println("\n→ Notice: Failures happen when generators are heavily loaded!")
end

println("\n" * "─"^70)
println("KEY INSIGHT: The coupling loop in action!")
println("High load → Heavy dispatch → High loading → Increased failures → More load shed")
println("─"^70)


# =============================================================================
# STEP 4: Understanding Rolling-Horizon vs Day-Ahead
# =============================================================================

println("\n\n### STEP 4: Rolling-Horizon UC - Why Look Ahead? ###\n")

"""
Compare different UC strategies:
1. Myopic (1-hour): Fast but may violate constraints
2. Day-ahead (24-hour): Good foresight but forecast error
3. Rolling-horizon (6-hour): Sweet spot
"""

println("Conceptual comparison:")
println()
println("┌─────────────────┬──────────┬───────────┬────────────────┐")
println("│ Strategy        │ Horizon  │ Accuracy  │ Computational  │")
println("├─────────────────┼──────────┼───────────┼────────────────┤")
println("│ Myopic (1-hour) │ t to t+1 │ Perfect   │ Very fast      │")
println("│                 │          │ boundary  │                │")
println("├─────────────────┼──────────┼───────────┼────────────────┤")
println("│ Rolling (6-hour)│ t to t+6 │ Good      │ Moderate       │")
println("│                 │          │ balance   │ (Paper choice) │")
println("├─────────────────┼──────────┼───────────┼────────────────┤")
println("│ Day-ahead (24h) │ t to t+24│ Forecast  │ Slow           │")
println("│                 │          │ error     │                │")
println("└─────────────────┴──────────┴───────────┴────────────────┘")
println()

println("Example: Why horizon matters")
println()
println("Suppose at 8 AM:")
println("  • Current load: 120 MW")
println("  • Generator ramp limit: 30 MW/hour")
println("  • Peak load at 2 PM (6 hours later): 200 MW")
println()

println("Myopic strategy (horizon=1):")
println("  → Only sees 120 MW now")
println("  → Commits cheap, low-capacity units")
println("  → At 2 PM: Can't ramp up fast enough → Load shed!")
println()

println("Rolling-horizon strategy (horizon=6):")
println("  → Sees the 2 PM peak coming")
println("  → Pre-commits larger units at 8 AM")
println("  → At 2 PM: Sufficient capacity → No load shed")
println()

println("Day-ahead strategy (horizon=24):")
println("  → Optimizes for entire day")
println("  → But: Uses forecast wind (uncertain)")
println("  → If actual wind < forecast → Load shed anyway")
println()

println("─"^70)
println("KEY INSIGHT: Rolling-horizon balances foresight and adaptation!")
println("Look far enough ahead to avoid commitment mistakes,")
println("but update frequently with actual renewable realizations.")
println("─"^70)


# =============================================================================
# STEP 5: Understanding Convergence - Why CV Matters
# =============================================================================

println("\n\n### STEP 5: Sequential MCS - Convergence Monitoring ###\n")

"""
Why do we need CV (Coefficient of Variation)?
Each scenario gives different results (randomness).
We need many scenarios until the MEAN stabilizes.
"""

println("Simulate convergence behavior:")
println()

# Generate fake LOLP samples (realistic distribution)
Random.seed!(42)
true_LOLP = 0.012  # Unknown true value

function generate_scenario_LOLP()
    # Each scenario samples failures randomly
    # LOLP will vary around true value
    return max(0.0, true_LOLP + 0.005 * randn())
end

LOLP_samples = [generate_scenario_LOLP() for _ in 1:200]

# Compute running statistics
function compute_running_stats(samples)
    n = length(samples)
    means = [mean(samples[1:i]) for i in 10:n]
    stds = [std(samples[1:i]) for i in 10:n]
    CVs = [std(samples[1:i]) / (mean(samples[1:i]) * sqrt(i)) for i in 10:n]
    return means, stds, CVs
end

means, stds, CVs = compute_running_stats(LOLP_samples)
scenarios = 10:length(LOLP_samples)

# Show key milestones
println("Scenario   Mean LOLP    Std Dev      CV        Status")
println("─"^60)
for i in [10, 25, 50, 100, 150, 200]
    idx = i - 9
    if idx <= length(means)
        cv = CVs[idx]
        status = cv <= 0.05 ? "✓ Converged" : "  Running..."
        println("$(lpad(i, 4))     $(round(means[idx], digits=4))     $(round(stds[idx], digits=4))     $(round(cv, digits=3))     $status")
    end
end

println()
println("Convergence achieved at scenario ~$(findfirst(x -> x <= 0.05, CVs) + 9)")

# Visualize convergence
plot(scenarios, CVs,
     xlabel="Number of Scenarios",
     ylabel="Coefficient of Variation",
     title="Convergence Monitoring",
     linewidth=2,
     label="CV",
     legend=:topright)
hline!([0.05], label="Target (5%)", linestyle=:dash, linewidth=2, color=:red)
savefig("step5_convergence.png")
println("\n✓ Plot saved as 'step5_convergence.png'")

println("\n" * "─"^70)
println("KEY INSIGHT: CV tells us when to stop simulating!")
println("CV = (Standard Error) / Mean")
println("CV < 0.05 means our estimate is within ~5% of true value (95% confidence)")
println("─"^70)


# =============================================================================
# SUMMARY
# =============================================================================

println("\n\n" * "="^70)
println("TUTORIAL COMPLETE - Summary")
println("="^70)
println()
println("You've learned the 5 key components:")
println()
println("1. DBBN: Non-parametric distributions capture tail events")
println("   → Better for reliability (we care about low-probability events)")
println()
println("2. PHM: Stress-dependent failure rates create realism")
println("   → Operational decisions affect reliability dynamically")
println()
println("3. Coupling: UC ↔ PHM creates bidirectional feedback")
println("   → Loading affects failures, failures affect next dispatch")
println()
println("4. Rolling-horizon: Balance foresight and adaptation")
println("   → Look ahead enough, but update with real data frequently")
println()
println("5. Convergence: CV monitors when estimate is accurate")
println("   → Stop simulating when CV < target (e.g., 5%)")
println()
println("="^70)
println()
println("Next steps:")
println("  • Review lecture slides (MachineLearning.tex) for mathematics")
println("  • Run full example: include(\"ML_Reliability_Examples.jl\")")
println("  • Read guide: ML_Reliability_Guide.md")
println()
println("Generated plots:")
println("  • step1_dbbn_samples.png - DBBN distribution")
println("  • step2_phm_loading.png - PHM stress effect")
println("  • step5_convergence.png - CV convergence")
println()
println("="^70)
