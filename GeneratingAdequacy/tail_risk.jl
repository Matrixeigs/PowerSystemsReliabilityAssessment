using Plots
using Printf
using Statistics
using Random

# Reuse existing structures
include("generating_adequacy_comprehensive.jl")

# ==============================================================================
# 1. DETAILED MONTE CARLO (Returns Arrays, not just Averages)
# ==============================================================================
function run_detailed_mc(gens::Vector{Generator}, base_load::Vector{Float64}, 
                        lfu_sigma_percent::Float64, n_years::Int)
    
    println("Running Detailed Monte Carlo ($n_years years)...")
    
    # Data Collection
    yearly_lole_distribution = Float64[] # LOLE for each of the 2000 years
    hourly_failure_prob = zeros(length(base_load)) # Probability of failure per hour
    
    total_hours = length(base_load)
    lfu_std_dev = maximum(base_load) * (lfu_sigma_percent / 100.0)
    gen_energy_state = zeros(Float64, length(gens))
    
    for year in 1:n_years
        fill!(gen_energy_state, 0.0)
        year_lole_count = 0.0
        
        for h in 1:total_hours
            current_week = div(h-1, 168) + 1
            
            # --- 1. Capacity State ---
            cap_unlimited = 0.0
            cap_elu_available = 0.0
            available_elu_indices = Int[]
            
            for (i, g) in enumerate(gens)
                # Maintenance
                if current_week >= g.scheduled_outage_start && 
                   current_week < g.scheduled_outage_start + g.maintenance_weeks
                    continue
                end
                # Random Failure
                if rand() < g.for_rate; continue; end
                # Energy Limit
                if g.energy_limit < Inf
                    if gen_energy_state[i] >= g.energy_limit
                        continue # Exhausted
                    else
                        cap_elu_available += g.capacity
                        push!(available_elu_indices, i)
                    end
                else
                    cap_unlimited += g.capacity
                end
            end
            
            # --- 2. Load & Dispatch ---
            actual_load = base_load[h] + randn() * lfu_std_dev
            unserved = max(0.0, actual_load - cap_unlimited)
            deficit = 0.0
            
            if unserved > 0
                if unserved > cap_elu_available
                    deficit = unserved - cap_elu_available
                    for idx in available_elu_indices
                        gen_energy_state[idx] += gens[idx].capacity
                    end
                else
                    needed = unserved
                    for idx in available_elu_indices
                        share = needed * (gens[idx].capacity / cap_elu_available)
                        gen_energy_state[idx] += share
                    end
                end
            end
            
            # --- 3. Record Data ---
            if deficit > 0
                year_lole_count += 1.0
                hourly_failure_prob[h] += 1.0
            end
        end
        push!(yearly_lole_distribution, year_lole_count)
    end
    
    # Normalize hourly counts to probabilities
    hourly_failure_prob ./= n_years
    
    return yearly_lole_distribution, hourly_failure_prob
end

# ==============================================================================
# 2. DETAILED ANALYTICAL (Extracts Hourly Risk Profile)
# ==============================================================================
function run_detailed_analytical(gens::Vector{Generator}, base_load::Vector{Float64}, 
                               lfu_sigma_percent::Float64)
    
    println("Running Detailed Analytical Profile...")
    hourly_risk_profile = zeros(length(base_load))
    
    step_size = 20.0
    lfu_mw = maximum(base_load) * (lfu_sigma_percent / 100.0)
    lfu_dist = get_lfu_distribution()
    
    # 1. Run Iterations to set effective_q
    for i in 1:5; update_elu!(gens, base_load, step_size, lfu_mw); end
    
    # 2. Calculate Hourly Risk
    for w in 1:52
        week_gens = Generator[]
        for g in gens
            if !(w >= g.scheduled_outage_start && w < g.scheduled_outage_start + g.maintenance_weeks)
                push!(week_gens, g)
            end
        end
        
        copt = COPT([0.0], [1.0])
        for g in week_gens; copt = add_unit(copt, g, step_size); end
        
        start_h = (w-1)*168 + 1
        end_h = min(w*168, 8760)
        installed = copt.capacity_outage[end]
        
        for h in start_h:end_h
            load = base_load[h]
            risk_h = 0.0
            for (z, p_z) in lfu_dist
                res = installed - (load + z*lfu_mw)
                for (i, out) in enumerate(copt.capacity_outage)
                    if out > res
                        risk_h += copt.probability[i] * p_z
                    end
                end
            end
            hourly_risk_profile[h] = risk_h
        end
    end
    
    return sum(hourly_risk_profile), hourly_risk_profile
end

# ==============================================================================
# 3. VISUALIZATION SCRIPT
# ==============================================================================
function visualize_tail_risk()
    # --- Setup ---
    gens = [
        Generator("Nuclear", 400.0, 0.02, 4, Inf),
        Generator("Coal_A",  300.0, 0.04, 3, Inf),
        Generator("Coal_B",  300.0, 0.04, 3, Inf),
        Generator("Gas",     150.0, 0.05, 2, Inf),
        Generator("Hydro_ELU", 200.0, 0.01, 2, 200.0 * 50.0), # 50 Hours Limit
        Generator("Old_56",   56.0, 0.10, 0, Inf)
    ]
    hours = 1:8760
    base_load = [750.0 + 300.0 * sin((h-2000)/8760 * 2π) + 50.0*randn() for h in hours]
    base_load = max.(0.0, base_load)
    weekly_peaks = [maximum(base_load[((w-1)*168+1):min(w*168,8760)]) for w in 1:52]
    schedule_maintenance!(gens, weekly_peaks)
    
    # --- Run Calculations ---
    ana_total, ana_profile = run_detailed_analytical(gens, base_load, 5.0)
    mc_dist, mc_profile = run_detailed_mc(gens, base_load, 5.0, 2000)
    
    # --- PLOT 1: The "Tail Risk" Histogram ---
    # This shows WHY the MC average is so high.
    p1 = histogram(mc_dist, bins=50, normalize=:pdf, alpha=0.6, color=:blue,
                   label="MC Yearly Outcomes", title="1. Distribution of Annual Risk (Tail Risk)",
                   xlabel="LOLE (Hours/Year)", ylabel="Probability Density")
    
    # Add Analytical Mean Line
    vline!([ana_total], color=:red, linewidth=3, label="Analytical Prediction")
    # Add MC Mean Line
    vline!([mean(mc_dist)], color=:blue, linestyle=:dash, linewidth=3, label="MC Average")
    
    # --- PLOT 2: The "Mechanism" (Zoomed Peak Week) ---
    # Find the week with the highest risk to zoom in
    peak_hour_idx = argmax(mc_profile)
    window = (peak_hour_idx - 50):(peak_hour_idx + 50)
    
    p2 = plot(window, mc_profile[window], fill=(0, 0.3, :blue), color=:blue,
              label="MC Risk (Sequential)", title="2. Hourly Risk Profile (Peak Window)",
              xlabel="Hour of Year", ylabel="Probability of Loss")
    
    plot!(window, ana_profile[window], color=:red, linewidth=2, 
          label="Analytical Risk (Smoothed)")
          
    # --- PLOT 3: Cumulative Risk (The Gap Accumulation) ---
    # This shows how the difference builds up over the year
    p3 = plot(cumsum(mc_profile), label="MC Cumulative LOLE", color=:blue, linewidth=2,
              title="3. Accumulation of Risk over Year", xlabel="Hour", ylabel="Cumulative LOLE")
    plot!(cumsum(ana_profile), label="Analytical Cumulative LOLE", color=:red, linestyle=:dash, linewidth=2)

    display(p1)
    display(p2)
    display(p3)
    
    println("\nVisualization Complete.")
    println("Plot 1: Notice if the Red Line (Analytical) is far to the left of the Blue Histogram.")
    println("Plot 2: Notice if MC has sharp spikes (Hydro Exhaustion) vs Analytical smooth curves.")
end

visualize_tail_risk()
