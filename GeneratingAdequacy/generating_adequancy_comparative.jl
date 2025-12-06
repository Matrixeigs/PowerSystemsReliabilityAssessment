using Plots
using Printf
using Statistics
using Random

# 1. REUSE EXISTING METHODS
# This imports Generator, COPT, add_unit, update_elu!, schedule_maintenance!
# Note: This might run the simulation in the original file once upon loading.
include("generating_adequacy_comprehensive.jl")

# ==============================================================================
# 2. NEW MONTE CARLO METHOD (Using External State)
# ==============================================================================

function run_monte_carlo_simulation(gens::Vector{Generator}, base_load::Vector{Float64}, 
                                  lfu_sigma_percent::Float64, n_years::Int)
    
    println("\n--- Starting Monte Carlo Simulation ($n_years years) ---")
    
    # METRICS
    total_lole_hours = 0.0
    total_eue_mwh = 0.0
    lole_history = Float64[]
    
    # CONSTANTS
    total_hours = length(base_load)
    lfu_std_dev = maximum(base_load) * (lfu_sigma_percent / 100.0)
    
    # EXTERNAL STATE ARRAYS
    # Since 'Generator' struct doesn't have 'current_energy', we track it here.
    # Map: gen_index -> energy_used_mwh
    gen_energy_state = zeros(Float64, length(gens))
    
    for year in 1:n_years
        # Reset Energy State for new year
        fill!(gen_energy_state, 0.0)
        
        year_lole = 0.0
        year_eue = 0.0
        
        for h in 1:total_hours
            current_week = div(h-1, 168) + 1
            
            # 1. Identify Available Capacity
            cap_unlimited = 0.0
            cap_elu_available = 0.0
            
            # We store indices of available ELUs to dispatch them later
            available_elu_indices = Int[]
            
            for (i, g) in enumerate(gens)
                # A. Maintenance Check (Reusing the field from original struct)
                if current_week >= g.scheduled_outage_start && 
                   current_week < g.scheduled_outage_start + g.maintenance_weeks
                    continue
                end
                
                # B. Random Failure (Monte Carlo Sampling)
                if rand() < g.for_rate
                    continue
                end
                
                # C. Energy Limit Check (Using External State)
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
            
            # 2. Load Uncertainty
            actual_load = base_load[h] + randn() * lfu_std_dev
            
            # 3. Dispatch Logic
            unserved_by_thermal = max(0.0, actual_load - cap_unlimited)
            deficit = 0.0
            
            if unserved_by_thermal > 0
                if unserved_by_thermal > cap_elu_available
                    # System Failure: Even with all ELUs, we can't meet load
                    deficit = unserved_by_thermal - cap_elu_available
                    
                    # All available ELUs run at full capacity
                    for idx in available_elu_indices
                        gen_energy_state[idx] += gens[idx].capacity
                    end
                else
                    # ELUs can cover it. Dispatch proportionally.
                    needed = unserved_by_thermal
                    for idx in available_elu_indices
                        g = gens[idx]
                        # Share proportional to capacity
                        share = needed * (g.capacity / cap_elu_available)
                        gen_energy_state[idx] += share
                    end
                end
            end
            
            # 4. Record Indices
            if deficit > 0
                year_lole += 1.0
                year_eue += deficit
            end
        end
        
        total_lole_hours += year_lole
        total_eue_mwh += year_eue
        
        if year % 100 == 0
            push!(lole_history, total_lole_hours / year)
        end
    end
    
    return total_lole_hours / n_years, total_eue_mwh / n_years, lole_history
end

# ==============================================================================
# 3. COMPARATIVE ANALYSIS DRIVER
# ==============================================================================

function perform_comparison()
    println("\n==================================================")
    println("   STARTING COMPARATIVE ANALYSIS")
    println("==================================================")

    # --- A. SETUP INPUTS (Identical to Original) ---
    gens = [
        Generator("Nuclear", 400.0, 0.02, 4, Inf),
        Generator("Coal_A",  300.0, 0.04, 3, Inf),
        Generator("Coal_B",  300.0, 0.04, 3, Inf),
        Generator("Gas",     150.0, 0.05, 2, Inf),
        # The Critical ELU: 200 MW with only 50 hours of water
        Generator("Hydro_ELU", 200.0, 0.01, 2, 200.0 * 50.0), 
        Generator("Old_56",   56.0, 0.10, 0, Inf)
    ]
    
    hours = 1:8760
    base_load = [750.0 + 300.0 * sin((h-2000)/8760 * 2π) + 50.0*randn() for h in hours]
    base_load = max.(0.0, base_load)
    weekly_peaks = [maximum(base_load[((w-1)*168+1):min(w*168,8760)]) for w in 1:52]
    lfu_mw = maximum(base_load) * 0.05 # 5% LFU

    # --- B. REUSE ANALYTICAL METHODS ---
    println("\n[1] Running Analytical Method (Reusing Code)...")
    
    # 1. Reuse Schedule Logic
    schedule_maintenance!(gens, weekly_peaks)
    
    # 2. Reuse ELU Update Logic (Iterative)
    # We run this to ensure the 'effective_q' is updated for the analytical calculation
    println("    Running ELU Iterations...")
    for i in 1:5
        update_elu!(gens, base_load, 20.0, lfu_mw)
    end
    
    # 3. Calculate Final Analytical LOLE (Manual Loop using reused 'add_unit')
    # We re-implement the risk summation loop here to get the exact number
    ana_lole = 0.0
    lfu_dist = get_lfu_distribution() # Reused function
    
    for w in 1:52
        # Filter maintenance
        week_gens = Generator[]
        for g in gens
            if !(w >= g.scheduled_outage_start && w < g.scheduled_outage_start + g.maintenance_weeks)
                push!(week_gens, g)
            end
        end
        
        # Build COPT (Reused function)
        copt = COPT([0.0], [1.0])
        for g in week_gens; copt = add_unit(copt, g, 20.0); end
        
        # Calc Risk
        start_h, end_h = (w-1)*168+1, min(w*168, 8760)
        installed = copt.capacity_outage[end]
        
        for h in start_h:end_h
            load = base_load[h]
            for (z, p_z) in lfu_dist
                res = installed - (load + z*lfu_mw)
                for (i, out) in enumerate(copt.capacity_outage)
                    if out > res
                        ana_lole += copt.probability[i] * p_z
                    end
                end
            end
        end
    end

    # --- C. RUN MONTE CARLO METHOD ---
    println("\n[2] Running Monte Carlo Method (New Code)...")
    # Note: MC uses the same 'gens', but ignores 'effective_q'. 
    # It uses 'energy_limit' + 'gen_energy_state' directly.
    mc_lole, mc_eue, history = run_monte_carlo_simulation(gens, base_load, 5.0, 2000)

    # --- D. REPORT ---
    println("\n--------------------------------------------------")
    println("FINAL RESULTS COMPARISON")
    println("--------------------------------------------------")
    @printf("Analytical LOLE (Iterative ELU): %.4f hours/year\n", ana_lole)
    @printf("Monte Carlo LOLE (Sequential):   %.4f hours/year\n", mc_lole)
    
    diff = abs(mc_lole - ana_lole)
    println("\nConclusion:")
    if diff < 50.0
        println("SUCCESS: The methods match closely! The Iterative Analytical method successfully approximated the complex ELU behavior.")
    else
        println("NOTICE: There is a gap. This highlights the 'Tail Risk' that Monte Carlo captures better than convolution.")
    end

    # Plot
    p = plot(100:100:2000, history, label="Monte Carlo", 
             xlabel="Years", ylabel="LOLE", title="Method Comparison", lw=2)
    hline!([ana_lole], label="Analytical (Iterative)", color=:red, linestyle=:dash, lw=2)
    display(p)
end

perform_comparison()
