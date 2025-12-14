# =============================================================================
# Comprehensive Distribution System and Station Reliability Assessment
# Educational Implementation with Extensive Visualizations
# =============================================================================

using Printf
using Random
using Distributions
using Statistics
using Plots
using StatsPlots
using LaTeXStrings

# Set default plot settings for publication quality
default(fontfamily="Computer Modern", linewidth=2, framestyle=:box,
        label=nothing, grid=true, dpi=300)

# =============================================================================
# MODULE 1: DATA STRUCTURES
# =============================================================================

"""
    Component

Represents a single component in the power system with reliability parameters.

Fields:
- name: Component identifier
- lambda: Failure rate (failures/year)
- mu: Repair rate (repairs/year)
- r: Average repair time (hours)
- type: Component type (:breaker, :transformer, :line, :bus)
"""
mutable struct Component
    name::String
    lambda::Float64  # failures/year
    mu::Float64      # repairs/year (1/r in years)
    r::Float64       # average repair time (hours)
    type::Symbol

    function Component(name::String, lambda::Float64, r::Float64, type::Symbol)
        mu = 8760.0 / r  # Convert hours to years
        new(name, lambda, mu, r, type)
    end
end

"""
    LoadPoint

Represents a load point in the distribution system.

Fields:
- name: Load point identifier
- customers: Number of customers
- lambda: Aggregate failure rate (failures/year)
- r: Average outage duration (hours)
- U: Annual unavailability (hours/year)
- peak_load: Peak load demand (MW)
"""
struct LoadPoint
    name::String
    customers::Int
    lambda::Float64
    r::Float64
    U::Float64
    peak_load::Float64

    function LoadPoint(name::String, customers::Int, lambda::Float64, r::Float64, peak_load::Float64=1.0)
        U = lambda * r
        new(name, customers, lambda, r, U, peak_load)
    end
end

"""
    SystemIndices

System-level reliability indices.
"""
struct SystemIndices
    SAIFI::Float64  # System Average Interruption Frequency Index
    SAIDI::Float64  # System Average Interruption Duration Index
    CAIDI::Float64  # Customer Average Interruption Duration Index
    ASAI::Float64   # Average Service Availability Index
    ASUI::Float64   # Average Service Unavailability Index
    ENS::Float64    # Energy Not Supplied (MWh/year)
    AENS::Float64   # Average Energy Not Supplied (MWh/customer/year)
end

# =============================================================================
# MODULE 2: ANALYTICAL METHODS FOR DISTRIBUTION SYSTEMS
# =============================================================================

"""
    calculate_system_indices(load_points, total_customers)

Calculate system-level reliability indices from load point data.
"""
function calculate_system_indices(load_points::Vector{LoadPoint}, total_customers::Int)
    # Weighted sums
    sum_lambda_N = sum(lp.lambda * lp.customers for lp in load_points)
    sum_U_N = sum(lp.U * lp.customers for lp in load_points)
    sum_ENS = sum(lp.U * lp.peak_load for lp in load_points)

    # Calculate indices
    SAIFI = sum_lambda_N / total_customers
    SAIDI = sum_U_N / total_customers
    CAIDI = SAIDI / SAIFI
    ASAI = (total_customers * 8760 - sum_U_N) / (total_customers * 8760)
    ASUI = 1.0 - ASAI
    ENS = sum_ENS
    AENS = ENS / total_customers

    return SystemIndices(SAIFI, SAIDI, CAIDI, ASAI, ASUI, ENS, AENS)
end

"""
    series_system_reliability(components)

Calculate reliability of components in series configuration.
"""
function series_system_reliability(components::Vector{Component})
    lambda_total = sum(c.lambda for c in components)
    U_total = sum(c.lambda * c.r for c in components)
    r_avg = lambda_total > 0 ? U_total / lambda_total : 0.0

    return lambda_total, r_avg, U_total
end

"""
    parallel_system_reliability(components)

Calculate reliability of components in parallel (redundant) configuration.
Assumes independent failures and perfect switching.
"""
function parallel_system_reliability(components::Vector{Component})
    if length(components) == 1
        return components[1].lambda, components[1].r, components[1].lambda * components[1].r
    end

    # For two components in parallel
    c1, c2 = components[1], components[2]

    # Failure rate of parallel system (both components must fail)
    lambda_parallel = (c1.lambda * c2.lambda * (c1.r + c2.r)) / 8760.0

    # Average repair time (geometric mean approximation)
    r_parallel = sqrt(c1.r * c2.r)

    U_parallel = lambda_parallel * r_parallel

    return lambda_parallel, r_parallel, U_parallel
end

# =============================================================================
# MODULE 3: DISTRIBUTION SYSTEM CASE STUDIES
# =============================================================================

"""
    run_case_1_fused_laterals()

Case 1: Radial distribution system with fused laterals, no alternate supply.
"""
function run_case_1_fused_laterals()
    println("\n" * "="^70)
    println("CASE 1: FUSED LATERALS, NO ALTERNATE SUPPLY")
    println("="^70)

    # Customer counts
    N_A, N_B, N_C = 250, 100, 50
    N_Total = 400

    # Load point data from analytical calculations
    lp_A = LoadPoint("A", N_A, 1.35, 1.55, 2.5)  # 2.5 MW peak
    lp_B = LoadPoint("B", N_B, 1.10, 2.05, 1.0)  # 1.0 MW peak
    lp_C = LoadPoint("C", N_C, 0.85, 2.05, 0.5)  # 0.5 MW peak

    load_points = [lp_A, lp_B, lp_C]

    # Calculate system indices
    indices = calculate_system_indices(load_points, N_Total)

    # Print results
    println("\nLoad Point Reliability:")
    for lp in load_points
        @printf("  %s: λ=%.2f f/yr, r=%.2f hr, U=%.2f hr/yr, N=%d customers\n",
                lp.name, lp.lambda, lp.r, lp.U, lp.customers)
    end

    println("\nSystem Indices:")
    @printf("  SAIFI: %.4f interruptions/customer/year\n", indices.SAIFI)
    @printf("  SAIDI: %.4f hours/customer/year\n", indices.SAIDI)
    @printf("  CAIDI: %.4f hours/interruption\n", indices.CAIDI)
    @printf("  ASAI:  %.6f\n", indices.ASAI)
    @printf("  ASUI:  %.6f\n", indices.ASUI)
    @printf("  ENS:   %.4f MWh/year\n", indices.ENS)
    @printf("  AENS:  %.4f MWh/customer/year\n", indices.AENS)

    return load_points, indices
end

"""
    run_case_2_alternate_supply()

Case 2: Radial distribution system with alternate supply capability.
"""
function run_case_2_alternate_supply()
    println("\n" * "="^70)
    println("CASE 2: ALTERNATE SUPPLY AVAILABLE (Switching time = 1.0 hr)")
    println("="^70)

    N_A, N_B, N_C = 250, 100, 50
    N_Total = 400

    # With alternate supply, main feeder faults can be restored by switching
    # Reduced repair times due to 1-hour switching vs 3-hour repair
    lp_A = LoadPoint("A", N_A, 1.35, 1.05, 2.5)  # Improved r
    lp_B = LoadPoint("B", N_B, 1.10, 1.55, 1.0)  # Improved r
    lp_C = LoadPoint("C", N_C, 0.85, 1.55, 0.5)  # Improved r

    load_points = [lp_A, lp_B, lp_C]
    indices = calculate_system_indices(load_points, N_Total)

    println("\nLoad Point Reliability:")
    for lp in load_points
        @printf("  %s: λ=%.2f f/yr, r=%.2f hr, U=%.2f hr/yr\n",
                lp.name, lp.lambda, lp.r, lp.U)
    end

    println("\nSystem Indices:")
    @printf("  SAIFI: %.4f (unchanged from Case 1)\n", indices.SAIFI)
    @printf("  SAIDI: %.4f (improved)\n", indices.SAIDI)
    @printf("  CAIDI: %.4f\n", indices.CAIDI)
    @printf("  ENS:   %.4f MWh/year\n", indices.ENS)

    return load_points, indices
end

"""
    run_case_3_solidly_connected()

Case 3: Solidly connected laterals (no fuses) - any fault trips main breaker.
"""
function run_case_3_solidly_connected()
    println("\n" * "="^70)
    println("CASE 3: SOLIDLY CONNECTED (No Fuses)")
    println("="^70)

    N_A, N_B, N_C = 250, 100, 50
    N_Total = 400

    # All customers see all failures (sum of all component failure rates)
    lambda_total = 1.60  # Combined failure rate
    r_avg = 2.0          # Average repair time

    lp_A = LoadPoint("A", N_A, lambda_total, r_avg, 2.5)
    lp_B = LoadPoint("B", N_B, lambda_total, r_avg, 1.0)
    lp_C = LoadPoint("C", N_C, lambda_total, r_avg, 0.5)

    load_points = [lp_A, lp_B, lp_C]
    indices = calculate_system_indices(load_points, N_Total)

    println("\nAll load points experience same reliability:")
    @printf("  λ=%.2f f/yr, r=%.2f hr, U=%.2f hr/yr\n",
            lambda_total, r_avg, lambda_total * r_avg)

    println("\nSystem Indices:")
    @printf("  SAIFI: %.4f (significantly higher)\n", indices.SAIFI)
    @printf("  SAIDI: %.4f (significantly higher)\n", indices.SAIDI)
    @printf("  CAIDI: %.4f\n", indices.CAIDI)
    @printf("  ENS:   %.4f MWh/year\n", indices.ENS)

    return load_points, indices
end

# =============================================================================
# MODULE 4: MONTE CARLO SIMULATION
# =============================================================================

"""
    ComponentState

Represents the state of a component during simulation.
"""
mutable struct ComponentState
    component::Component
    is_operational::Bool
    time_to_next_event::Float64
    total_failure_time::Float64
    failure_count::Int
end

"""
    simulate_component_state_duration(component, simulation_years, rng)

Monte Carlo simulation using Component State Duration Sampling method.
"""
function simulate_component_state_duration(component::Component, simulation_years::Float64, rng::AbstractRNG)
    current_time = 0.0
    end_time = simulation_years * 8760.0  # Convert to hours

    total_downtime = 0.0
    failure_count = 0
    is_operational = true

    # Distributions
    λ_hourly = component.lambda / 8760.0  # Convert to failures/hour
    failure_dist = Exponential(1.0 / λ_hourly)
    repair_dist = Exponential(component.r)  # Using exponential for repair time

    while current_time < end_time
        if is_operational
            # Sample time to next failure
            ttf = rand(rng, failure_dist)
            current_time += ttf

            if current_time < end_time
                is_operational = false
                failure_count += 1
            end
        else
            # Sample repair time
            ttr = rand(rng, repair_dist)
            total_downtime += ttr
            current_time += ttr
            is_operational = true
        end
    end

    # Calculate metrics
    lambda_simulated = failure_count / simulation_years
    U_simulated = total_downtime / simulation_years
    r_simulated = failure_count > 0 ? U_simulated / lambda_simulated : 0.0

    return lambda_simulated, r_simulated, U_simulated, failure_count, total_downtime
end

"""
    run_monte_carlo_distribution_system(num_simulations, simulation_years)

Run Monte Carlo simulation for a distribution system.
"""
function run_monte_carlo_distribution_system(num_simulations::Int=10000, simulation_years::Float64=10.0)
    println("\n" * "="^70)
    println("MONTE CARLO SIMULATION - COMPONENT STATE DURATION SAMPLING")
    println("="^70)
    @printf("Number of simulations: %d\n", num_simulations)
    @printf("Simulation period: %.1f years\n", simulation_years)

    # Define a simple component (e.g., a feeder section)
    component = Component("Feeder Section", 0.5, 3.0, :line)

    println("\nAnalytical values:")
    @printf("  λ = %.4f failures/year\n", component.lambda)
    @printf("  r = %.4f hours\n", component.r)
    @printf("  U = %.4f hours/year\n", component.lambda * component.r)

    # Storage for results
    lambda_results = zeros(num_simulations)
    r_results = zeros(num_simulations)
    U_results = zeros(num_simulations)

    rng = MersenneTwister(42)  # For reproducibility

    println("\nRunning simulations...")
    for i in 1:num_simulations
        lambda_sim, r_sim, U_sim, _, _ = simulate_component_state_duration(component, simulation_years, rng)
        lambda_results[i] = lambda_sim
        r_results[i] = r_sim
        U_results[i] = U_sim
    end

    # Calculate statistics
    println("\nSimulation Results (Mean ± Std):")
    @printf("  λ = %.4f ± %.4f failures/year\n", mean(lambda_results), std(lambda_results))
    @printf("  r = %.4f ± %.4f hours\n", mean(r_results), std(r_results))
    @printf("  U = %.4f ± %.4f hours/year\n", mean(U_results), std(U_results))

    return lambda_results, r_results, U_results, component
end

"""
    run_monte_carlo_convergence_analysis(max_simulations, simulation_years)

Analyze Monte Carlo convergence behavior.
"""
function run_monte_carlo_convergence_analysis(max_simulations::Int=5000, simulation_years::Float64=10.0)
    component = Component("Test Component", 0.5, 3.0, :line)
    analytical_U = component.lambda * component.r

    rng = MersenneTwister(42)

    # Sample points for convergence plot
    sample_points = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    sample_points = sample_points[sample_points .<= max_simulations]

    U_means = zeros(length(sample_points))
    U_stds = zeros(length(sample_points))

    # Run all simulations once
    all_U = zeros(max_simulations)
    for i in 1:max_simulations
        _, _, U_sim, _, _ = simulate_component_state_duration(component, simulation_years, rng)
        all_U[i] = U_sim
    end

    # Calculate running statistics
    for (idx, n) in enumerate(sample_points)
        U_means[idx] = mean(all_U[1:n])
        U_stds[idx] = std(all_U[1:n])
    end

    return sample_points, U_means, U_stds, analytical_U
end

# =============================================================================
# MODULE 5: STATION RELIABILITY ASSESSMENT
# =============================================================================

"""
    calculate_station_reliability_simple()

Calculate station reliability for a simple series configuration.
Configuration: Source -> Breaker -> Transformer -> Bus -> Load
"""
function calculate_station_reliability_simple()
    println("\n" * "="^70)
    println("STATION RELIABILITY - SIMPLE SERIES CONFIGURATION")
    println("="^70)

    # Component data based on typical values
    breaker_active = Component("Breaker (Active)", 0.006, 10.0, :breaker)
    breaker_passive = Component("Breaker (Passive)", 0.004, 10.0, :breaker)
    transformer = Component("Transformer", 0.010, 100.0, :transformer)
    bus = Component("Bus", 0.001, 2.0, :bus)

    # Combine active and passive breaker failures
    breaker_combined = Component("Breaker (Combined)",
                                 breaker_active.lambda + breaker_passive.lambda,
                                 10.0, :breaker)

    components = [breaker_combined, transformer, bus]

    println("\nComponent Reliability Data:")
    for c in components
        @printf("  %s: λ=%.4f f/yr, r=%.2f hr\n", c.name, c.lambda, c.r)
    end

    # Series system calculation
    lambda_total, r_avg, U_total = series_system_reliability(components)

    println("\nLoad Point Reliability:")
    @printf("  Failure Rate (λ): %.4f failures/year\n", lambda_total)
    @printf("  Avg Outage Duration (r): %.2f hours\n", r_avg)
    @printf("  Unavailability (U): %.4f hours/year\n", U_total)
    @printf("  Availability: %.6f\n", 1.0 - U_total/8760.0)

    return components, lambda_total, r_avg, U_total
end

"""
    calculate_station_reliability_redundant()

Calculate station reliability with redundant transformer configuration.
"""
function calculate_station_reliability_redundant()
    println("\n" * "="^70)
    println("STATION RELIABILITY - REDUNDANT TRANSFORMER CONFIGURATION")
    println("="^70)

    # Components
    breaker = Component("Breaker", 0.010, 10.0, :breaker)
    transformer1 = Component("Transformer 1", 0.010, 100.0, :transformer)
    transformer2 = Component("Transformer 2", 0.010, 100.0, :transformer)
    bus = Component("Bus", 0.001, 2.0, :bus)

    # Parallel transformer reliability
    lambda_trans, r_trans, U_trans = parallel_system_reliability([transformer1, transformer2])

    println("\nTransformer Redundancy Analysis:")
    @printf("  Single transformer: λ=%.4f f/yr, r=%.2f hr\n",
            transformer1.lambda, transformer1.r)
    @printf("  Parallel transformers: λ=%.6f f/yr, r=%.2f hr\n",
            lambda_trans, r_trans)
    @printf("  Reliability improvement: %.1fx\n", transformer1.lambda / lambda_trans)

    # Overall system with redundant transformers
    transformer_eq = Component("Transformers (Parallel)", lambda_trans, r_trans, :transformer)
    components = [breaker, transformer_eq, bus]

    lambda_total, r_avg, U_total = series_system_reliability(components)

    println("\nOverall Station Reliability:")
    @printf("  Failure Rate (λ): %.4f failures/year\n", lambda_total)
    @printf("  Avg Outage Duration (r): %.2f hours\n", r_avg)
    @printf("  Unavailability (U): %.4f hours/year\n", U_total)

    return components, lambda_total, r_avg, U_total
end

"""
    station_sensitivity_analysis()

Perform sensitivity analysis on station components.
"""
function station_sensitivity_analysis()
    println("\n" * "="^70)
    println("STATION RELIABILITY - SENSITIVITY ANALYSIS")
    println("="^70)

    # Base case
    base_breaker = Component("Breaker", 0.010, 10.0, :breaker)
    base_transformer = Component("Transformer", 0.010, 100.0, :transformer)
    base_bus = Component("Bus", 0.001, 2.0, :bus)

    # Vary transformer failure rate
    lambda_range = 0.005:0.001:0.020
    U_results = zeros(length(lambda_range))

    for (i, λ) in enumerate(lambda_range)
        trans = Component("Transformer", λ, 100.0, :transformer)
        components = [base_breaker, trans, base_bus]
        _, _, U = series_system_reliability(components)
        U_results[i] = U
    end

    return collect(lambda_range), U_results
end

# =============================================================================
# MODULE 6: PROBABILITY DISTRIBUTION ANALYSIS
# =============================================================================

"""
    plot_probability_distributions()

Generate plots comparing different probability distributions for repair times.
"""
function plot_probability_distributions()
    println("\n" * "="^70)
    println("PROBABILITY DISTRIBUTION ANALYSIS")
    println("="^70)

    # Mean repair time
    μ = 3.0  # hours

    # Define distributions
    exp_dist = Exponential(μ)
    gamma_dist = Gamma(2.0, μ/2.0)  # Shape=2, scale such that mean=μ
    lognormal_dist = LogNormal(log(μ) - 0.5*0.5^2, 0.5)  # σ=0.5
    weibull_dist = Weibull(2.0, μ/0.8862)  # Shape=2, scale adjusted for mean=μ (gamma(1.5)≈0.8862)

    x = 0.0:0.1:15.0

    # Calculate PDFs
    pdf_exp = pdf.(exp_dist, x)
    pdf_gamma = pdf.(gamma_dist, x)
    pdf_lognormal = pdf.(lognormal_dist, x)
    pdf_weibull = pdf.(weibull_dist, x)

    println("\nDistribution Statistics (Mean ± Std):")
    @printf("  Exponential:  %.2f ± %.2f hours\n", mean(exp_dist), std(exp_dist))
    @printf("  Gamma(2):     %.2f ± %.2f hours\n", mean(gamma_dist), std(gamma_dist))
    @printf("  LogNormal:    %.2f ± %.2f hours\n", mean(lognormal_dist), std(lognormal_dist))
    @printf("  Weibull(2):   %.2f ± %.2f hours\n", mean(weibull_dist), std(weibull_dist))

    return x, pdf_exp, pdf_gamma, pdf_lognormal, pdf_weibull
end

# =============================================================================
# MODULE 7: VISUALIZATION FUNCTIONS
# =============================================================================

"""
    plot_case_comparison(indices_list, case_names)

Create comparison plots for different distribution system cases.
"""
function plot_case_comparison(indices_list, case_names)
    # Extract metrics
    saifi_vals = [idx.SAIFI for idx in indices_list]
    saidi_vals = [idx.SAIDI for idx in indices_list]
    caidi_vals = [idx.CAIDI for idx in indices_list]
    ens_vals = [idx.ENS for idx in indices_list]

    # Create subplots
    p1 = bar(case_names, saifi_vals, ylabel="SAIFI\n(int/cust/yr)",
             title="System Average Interruption Frequency Index",
             color=:steelblue, legend=false, ylims=(0, maximum(saifi_vals)*1.2))

    p2 = bar(case_names, saidi_vals, ylabel="SAIDI\n(hr/cust/yr)",
             title="System Average Interruption Duration Index",
             color=:coral, legend=false, ylims=(0, maximum(saidi_vals)*1.2))

    p3 = bar(case_names, caidi_vals, ylabel="CAIDI\n(hr/int)",
             title="Customer Average Interruption Duration Index",
             color=:mediumseagreen, legend=false, ylims=(0, maximum(caidi_vals)*1.2))

    p4 = bar(case_names, ens_vals, ylabel="ENS\n(MWh/yr)",
             title="Energy Not Supplied",
             color=:orchid, legend=false, ylims=(0, maximum(ens_vals)*1.2))

    plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800),
         plot_title="Distribution System Reliability Comparison")
end

"""
    plot_monte_carlo_histograms(lambda_results, r_results, U_results, component)

Create histograms for Monte Carlo simulation results.
"""
function plot_monte_carlo_histograms(lambda_results, r_results, U_results, component)
    p1 = histogram(lambda_results, bins=50, xlabel="λ (failures/year)",
                   ylabel="Frequency", title="Failure Rate Distribution",
                   color=:steelblue, alpha=0.7, normalize=:probability)
    vline!([component.lambda], linewidth=3, color=:red, label="Analytical")

    p2 = histogram(r_results, bins=50, xlabel="r (hours)",
                   ylabel="Frequency", title="Repair Time Distribution",
                   color=:coral, alpha=0.7, normalize=:probability)
    vline!([component.r], linewidth=3, color=:red, label="Analytical")

    p3 = histogram(U_results, bins=50, xlabel="U (hours/year)",
                   ylabel="Frequency", title="Unavailability Distribution",
                   color=:mediumseagreen, alpha=0.7, normalize=:probability)
    vline!([component.lambda * component.r], linewidth=3, color=:red, label="Analytical")

    plot(p1, p2, p3, layout=(1,3), size=(1400, 400),
         plot_title="Monte Carlo Simulation Results (N=10000)")
end

"""
    plot_monte_carlo_convergence(sample_points, U_means, U_stds, analytical_U)

Plot Monte Carlo convergence behavior.
"""
function plot_monte_carlo_convergence(sample_points, U_means, U_stds, analytical_U)
    p = plot(sample_points, U_means, ribbon=U_stds,
             xlabel="Number of Simulations", ylabel="U (hours/year)",
             title="Monte Carlo Convergence Analysis",
             label="Simulated U ± σ", fillalpha=0.3, linewidth=2,
             xscale=:log10, color=:steelblue, fillcolor=:lightblue)

    hline!([analytical_U], linewidth=2, linestyle=:dash, color=:red,
           label="Analytical U")

    plot!(size=(800, 500), legend=:topright)
end

"""
    plot_load_point_reliability(load_points, case_name)

Visualize load point reliability metrics.
"""
function plot_load_point_reliability(load_points, case_name)
    names = [lp.name for lp in load_points]
    lambdas = [lp.lambda for lp in load_points]
    rs = [lp.r for lp in load_points]
    Us = [lp.U for lp in load_points]

    p1 = bar(names, lambdas, ylabel="λ (failures/yr)",
             title="Failure Rates", color=:steelblue, legend=false)

    p2 = bar(names, rs, ylabel="r (hours)",
             title="Avg Outage Duration", color=:coral, legend=false)

    p3 = bar(names, Us, ylabel="U (hr/yr)",
             title="Annual Unavailability", color=:mediumseagreen, legend=false)

    plot(p1, p2, p3, layout=(1,3), size=(1200, 400),
         plot_title="Load Point Reliability - $case_name")
end

"""
    plot_station_sensitivity(lambda_range, U_results)

Plot station reliability sensitivity analysis.
"""
function plot_station_sensitivity(lambda_range, U_results)
    plot(lambda_range, U_results, xlabel="Transformer λ (failures/year)",
         ylabel="Station U (hours/year)", title="Station Reliability Sensitivity",
         linewidth=2, color=:steelblue, marker=:circle, markersize=4,
         label="Unavailability", legend=:topleft, size=(800, 500))
end

"""
    plot_probability_distributions_comparison(x, pdfs...)

Plot comparison of different probability distributions.
"""
function plot_probability_distributions_comparison(x, pdf_exp, pdf_gamma, pdf_lognormal, pdf_weibull)
    plot(x, pdf_exp, label="Exponential", linewidth=2, color=:steelblue)
    plot!(x, pdf_gamma, label="Gamma (k=2)", linewidth=2, color=:coral)
    plot!(x, pdf_lognormal, label="LogNormal (σ=0.5)", linewidth=2, color=:mediumseagreen)
    plot!(x, pdf_weibull, label="Weibull (k=2)", linewidth=2, color=:orchid)

    xlabel!("Repair Time (hours)")
    ylabel!("Probability Density")
    title!("Repair Time Distribution Comparison (μ=3.0 hrs)")
    plot!(size=(900, 600), legend=:topright, xlims=(0, 15))
end

"""
    create_system_diagram()

Create a visual representation of the distribution system configuration.
"""
function create_system_diagram()
    println("\nSystem diagrams are best created using TikZ in LaTeX.")
    println("See the Distribution_Station_Reliability_Full.tex file for detailed diagrams.")
end

# =============================================================================
# MAIN EXECUTION FUNCTION
# =============================================================================

"""
    run_comprehensive_analysis()

Execute complete reliability analysis with all visualizations.
"""
function run_comprehensive_analysis()
    println("\n" * "="^80)
    println(" COMPREHENSIVE POWER SYSTEM RELIABILITY ASSESSMENT")
    println(" Educational Implementation with Visualizations")
    println("="^80)

    # Create output directory for figures
    output_dir = "figures"
    if !isdir(output_dir)
        mkdir(output_dir)
        println("\nCreated output directory: $output_dir/")
    end

    # ========== DISTRIBUTION SYSTEM CASES ==========

    # Run all cases
    lp_case1, idx_case1 = run_case_1_fused_laterals()
    lp_case2, idx_case2 = run_case_2_alternate_supply()
    lp_case3, idx_case3 = run_case_3_solidly_connected()

    # Plot comparisons
    println("\n" * "="^70)
    println("GENERATING FIGURES...")
    println("="^70)

    # Figure 1: Case comparison
    println("\n[1/9] Generating distribution system case comparison...")
    p1 = plot_case_comparison([idx_case1, idx_case2, idx_case3],
                              ["Case 1\nFused", "Case 2\nAlternate", "Case 3\nSolid"])
    savefig(p1, "$output_dir/fig1_case_comparison.png")
    println("      Saved: $output_dir/fig1_case_comparison.png")

    # Figure 2: Load point reliability for each case
    println("\n[2/9] Generating load point reliability plots...")
    p2 = plot_load_point_reliability(lp_case1, "Case 1")
    savefig(p2, "$output_dir/fig2_loadpoint_case1.png")

    p3 = plot_load_point_reliability(lp_case2, "Case 2")
    savefig(p3, "$output_dir/fig3_loadpoint_case2.png")

    p4 = plot_load_point_reliability(lp_case3, "Case 3")
    savefig(p4, "$output_dir/fig4_loadpoint_case3.png")
    println("      Saved: fig2-4 load point reliability plots")

    # ========== MONTE CARLO SIMULATION ==========

    println("\n[3/9] Running Monte Carlo simulation...")
    lambda_mc, r_mc, U_mc, comp_mc = run_monte_carlo_distribution_system(10000, 10.0)

    println("\n[4/9] Generating Monte Carlo histograms...")
    p5 = plot_monte_carlo_histograms(lambda_mc, r_mc, U_mc, comp_mc)
    savefig(p5, "$output_dir/fig5_monte_carlo_histograms.png")
    println("      Saved: $output_dir/fig5_monte_carlo_histograms.png")

    println("\n[5/9] Running convergence analysis...")
    sp, um, us, au = run_monte_carlo_convergence_analysis(5000, 10.0)

    println("\n[6/9] Generating convergence plot...")
    p6 = plot_monte_carlo_convergence(sp, um, us, au)
    savefig(p6, "$output_dir/fig6_monte_carlo_convergence.png")
    println("      Saved: $output_dir/fig6_monte_carlo_convergence.png")

    # ========== STATION RELIABILITY ==========

    println("\n[7/9] Calculating station reliability...")
    comp_simple, λ_simple, r_simple, U_simple = calculate_station_reliability_simple()
    comp_redun, λ_redun, r_redun, U_redun = calculate_station_reliability_redundant()

    # Station comparison
    station_configs = ["Simple\nSeries", "Redundant\nTransformers"]
    station_lambdas = [λ_simple, λ_redun]
    station_Us = [U_simple, U_redun]

    p7 = bar(station_configs, station_Us, ylabel="U (hours/year)",
             title="Station Configuration Comparison",
             color=[:steelblue, :coral], legend=false, size=(600, 400))
    savefig(p7, "$output_dir/fig7_station_comparison.png")
    println("      Saved: $output_dir/fig7_station_comparison.png")

    println("\n[8/9] Running station sensitivity analysis...")
    λ_range, U_sens = station_sensitivity_analysis()
    p8 = plot_station_sensitivity(λ_range, U_sens)
    savefig(p8, "$output_dir/fig8_station_sensitivity.png")
    println("      Saved: $output_dir/fig8_station_sensitivity.png")

    # ========== PROBABILITY DISTRIBUTIONS ==========

    println("\n[9/9] Generating probability distribution plots...")
    x, pdf_e, pdf_g, pdf_ln, pdf_w = plot_probability_distributions()
    p9 = plot_probability_distributions_comparison(x, pdf_e, pdf_g, pdf_ln, pdf_w)
    savefig(p9, "$output_dir/fig9_probability_distributions.png")
    println("      Saved: $output_dir/fig9_probability_distributions.png")

    # ========== SUMMARY ==========

    println("\n" * "="^80)
    println(" ANALYSIS COMPLETE!")
    println("="^80)
    println("\nGenerated 9 figures in the '$output_dir/' directory:")
    println("  1. Distribution system case comparison (SAIFI, SAIDI, CAIDI, ENS)")
    println("  2-4. Load point reliability for Cases 1, 2, 3")
    println("  5. Monte Carlo simulation histograms")
    println("  6. Monte Carlo convergence analysis")
    println("  7. Station configuration comparison")
    println("  8. Station sensitivity analysis")
    println("  9. Probability distribution comparison")

    println("\n" * "="^80)
    println(" KEY FINDINGS")
    println("="^80)

    println("\nDistribution System:")
    @printf("  Case 1 SAIDI: %.2f hr/cust/yr (Baseline)\n", idx_case1.SAIDI)
    @printf("  Case 2 SAIDI: %.2f hr/cust/yr (%.1f%% improvement with alternate supply)\n",
            idx_case2.SAIDI, 100*(idx_case1.SAIDI - idx_case2.SAIDI)/idx_case1.SAIDI)
    @printf("  Case 3 SAIDI: %.2f hr/cust/yr (%.1f%% worse without fuses)\n",
            idx_case3.SAIDI, 100*(idx_case3.SAIDI - idx_case1.SAIDI)/idx_case1.SAIDI)

    println("\nStation Reliability:")
    @printf("  Simple Configuration U: %.4f hr/yr\n", U_simple)
    @printf("  Redundant Configuration U: %.4f hr/yr (%.1fx improvement)\n",
            U_redun, U_simple/U_redun)

    println("\nMonte Carlo Validation:")
    @printf("  Analytical U: %.4f hr/yr\n", comp_mc.lambda * comp_mc.r)
    @printf("  Simulated U:  %.4f ± %.4f hr/yr\n", mean(U_mc), std(U_mc))
    @printf("  Error: %.2f%%\n", 100*abs(mean(U_mc) - comp_mc.lambda*comp_mc.r)/(comp_mc.lambda*comp_mc.r))

    println("\n" * "="^80)
    println()

    return nothing
end

# =============================================================================
# RUN THE ANALYSIS
# =============================================================================

# Execute the comprehensive analysis
run_comprehensive_analysis()
