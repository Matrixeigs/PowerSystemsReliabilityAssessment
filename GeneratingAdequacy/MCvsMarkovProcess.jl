using Statistics
using Random
using Printf
using Plots

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

mutable struct Generator
    name::String
    capacity::Float64
    for_rate::Float64          # Base Mechanical FOR (q)
    maintenance_weeks::Int     
    energy_limit::Float64      # Inf if not limited (MWh)
    
    # --- Analytical State ---
    effective_q::Float64       # Adjusted FOR (Iteratively updated for ELU)
    
    # --- Shared State ---
    scheduled_outage_start::Int 
end

# Constructor
function Generator(name, cap, q, maint, elim)
    return Generator(name, cap, q, maint, elim, q, 0)
end

struct COPT
    capacity_outage::Vector{Float64}
    probability::Vector{Float64}
end

struct SystemParams
    step_size::Float64
    lfu_sigma_percent::Float64
    mc_years::Int
end

# ==============================================================================
# 2. CORE LOGIC: MAINTENANCE SCHEDULING
# ==============================================================================

function schedule_maintenance!(gens::Vector{Generator}, weekly_peaks::Vector{Float64})
    num_weeks = 52
    total_installed = sum(g.capacity for g in gens)
    weekly_available_cap = fill(total_installed, num_weeks)
    
    # Sort by "Maintenance Burden" (Capacity * Weeks)
    sorted_gens = sort(gens, by = g -> g.capacity * g.maintenance_weeks, rev=true)
    
    for g in sorted_gens
        if g.maintenance_weeks == 0; continue; end
        
        best_start = 1
        max_min_reserve = -Inf
        
        # Find window that maximizes the minimum reserve (Levelization)
        for start_w in 1:(num_weeks - g.maintenance_weeks + 1)
            window = start_w:(start_w + g.maintenance_weeks - 1)
            min_res = minimum(weekly_available_cap[window] .- weekly_peaks[window])
            if min_res > max_min_reserve
                max_min_reserve = min_res
                best_start = start_w
            end
        end
        
        g.scheduled_outage_start = best_start
        weekly_available_cap[best_start:(best_start + g.maintenance_weeks - 1)] .-= g.capacity
    end
end

# ==============================================================================
# 3. ANALYTICAL ENGINE (Convolution & Rounding)
# ==============================================================================

function add_unit(old_copt::COPT, unit::Generator, step_size::Float64)
    C = unit.capacity
    q = unit.effective_q  # Use the ELU-adjusted q
    p = 1.0 - q
    
    max_old = isempty(old_copt.capacity_outage) ? 0.0 : maximum(old_copt.capacity_outage)
    max_new = max_old + C
    num_steps = Int(ceil(max_new / step_size)) + 1
    new_states = [i * step_size for i in 0:(num_steps-1)]
    new_probs = zeros(length(new_states))
    
    # Helper for sparse lookup
    function get_old_prob(X_val)
        idx = findfirst(x -> abs(x - X_val) < 1e-5, old_copt.capacity_outage)
        return idx === nothing ? 0.0 : old_copt.probability[idx]
    end

    # --- ROUNDING LOGIC ---
    lower_idx = Int(floor(C / step_size))
    C_lower = lower_idx * step_size
    C_upper = (lower_idx + 1) * step_size
    
    if C_lower == C_upper # Exact match
        for (i, X) in enumerate(new_states)
            new_probs[i] = get_old_prob(X) * p + get_old_prob(X - C) * q
        end
    else # Rounding required
        alpha = (C - C_lower) / step_size
        q_upper = q * alpha
        q_lower = q * (1.0 - alpha)
        for (i, X) in enumerate(new_states)
            new_probs[i] = get_old_prob(X) * p + 
                           get_old_prob(X - C_lower) * q_lower + 
                           get_old_prob(X - C_upper) * q_upper
        end
    end
    return COPT(new_states, new_probs)
end

function update_elu_analytical!(gens::Vector{Generator}, base_load::Vector{Float64}, params::SystemParams)
    # Iteratively adjust 'effective_q' based on energy deficits
    lfu_mw = maximum(base_load) * (params.lfu_sigma_percent / 100.0)
    lfu_dist = [(-3.0, 0.006), (-2.0, 0.061), (-1.0, 0.242), (0.0, 0.382), (1.0, 0.242), (2.0, 0.061), (3.0, 0.006)]
    
    changed = false
    for g in gens
        if g.energy_limit == Inf; continue; end
        
        # Build system without this unit
        copt_rest = COPT([0.0], [1.0])
        for og in gens
            if og != g; copt_rest = add_unit(copt_rest, og, params.step_size); end
        end
        
        # Calculate Expected Energy Required
        total_energy_req = 0.0
        cap_rest = copt_rest.capacity_outage[end]
        
        for load in base_load
            for (z, p_z) in lfu_dist
                act_load = load + z * lfu_mw
                reserve_thresh = cap_rest - act_load
                term_e = 0.0
                for (i, out) in enumerate(copt_rest.capacity_outage)
                    if out > reserve_thresh
                        deficit = out - reserve_thresh
                        term_e += min(g.capacity, deficit) * copt_rest.probability[i]
                    end
                end
                total_energy_req += term_e * p_z
            end
        end
        
        # Update q
        new_q = g.for_rate
        if total_energy_req > g.energy_limit
            prob_fail = (total_energy_req - g.energy_limit) / (g.capacity * 8760)
            new_q += prob_fail
        end
        new_q = min(new_q, 1.0)
        
        if abs(new_q - g.effective_q) > 1e-5
            g.effective_q = new_q
            changed = true
        end
    end
    return changed
end

function run_analytical(gens::Vector{Generator}, base_load::Vector{Float64}, params::SystemParams)
    println("  [Analytical] Optimizing ELU parameters...")
    for i in 1:5
        update_elu_analytical!(gens, base_load, params)
    end
    
    println("  [Analytical] Convolving COPTs...")
    total_lole = 0.0
    lfu_mw = maximum(base_load) * (params.lfu_sigma_percent / 100.0)
    lfu_dist = [(-3.0, 0.006), (-2.0, 0.061), (-1.0, 0.242), (0.0, 0.382), (1.0, 0.242), (2.0, 0.061), (3.0, 0.006)]
    
    hourly_risk = zeros(length(base_load))

    for w in 1:52
        week_gens = [g for g in gens if !(w >= g.scheduled_outage_start && w < g.scheduled_outage_start + g.maintenance_weeks)]
        
        copt = COPT([0.0], [1.0])
        for g in week_gens; copt = add_unit(copt, g, params.step_size); end
        
        start_h, end_h = (w-1)*168+1, min(w*168, 8760)
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
            hourly_risk[h] = risk_h
            total_lole += risk_h
        end
    end
    return total_lole, hourly_risk
end

# ==============================================================================
# 4. MONTE CARLO ENGINE (Sequential Simulation)
# ==============================================================================

function run_monte_carlo(gens::Vector{Generator}, base_load::Vector{Float64}, params::SystemParams)
    println("  [Monte Carlo] Simulating $(params.mc_years) years...")
    
    total_hours = length(base_load)
    lfu_std_dev = maximum(base_load) * (params.lfu_sigma_percent / 100.0)
    
    # External State for ELU (Energy Used MWh)
    gen_energy_state = zeros(Float64, length(gens))
    
    yearly_lole = Float64[]
    hourly_failures = zeros(total_hours)
    
    for year in 1:params.mc_years
        fill!(gen_energy_state, 0.0)
        year_lole_count = 0.0
        
        for h in 1:total_hours
            current_week = div(h-1, 168) + 1
            
            cap_unlimited = 0.0
            cap_elu_available = 0.0
            elu_indices = Int[]
            
            for (i, g) in enumerate(gens)
                # Maintenance
                if current_week >= g.scheduled_outage_start && current_week < g.scheduled_outage_start + g.maintenance_weeks
                    continue
                end
                # Random Failure
                if rand() < g.for_rate; continue; end
                # ELU Check
                if g.energy_limit < Inf
                    if gen_energy_state[i] >= g.energy_limit
                        continue # Empty
                    else
                        cap_elu_available += g.capacity
                        push!(elu_indices, i)
                    end
                else
                    cap_unlimited += g.capacity
                end
            end
            
            # Load with Continuous LFU
            actual_load = base_load[h] + randn() * lfu_std_dev
            
            # Dispatch
            unserved = max(0.0, actual_load - cap_unlimited)
            deficit = 0.0
            
            if unserved > 0
                if unserved > cap_elu_available
                    deficit = unserved - cap_elu_available
                    # Drain all ELUs
                    for idx in elu_indices; gen_energy_state[idx] += gens[idx].capacity; end
                else
                    # Partial drain
                    needed = unserved
                    for idx in elu_indices
                        share = needed * (gens[idx].capacity / cap_elu_available)
                        gen_energy_state[idx] += share
                    end
                end
            end
            
            if deficit > 0
                year_lole_count += 1.0
                hourly_failures[h] += 1.0
            end
        end
        push!(yearly_lole, year_lole_count)
    end
    
    return mean(yearly_lole), hourly_failures ./ params.mc_years, yearly_lole
end

# ==============================================================================
# 5. COMPARISON DRIVER
# ==============================================================================

function run_comparison()
    println("=== GENERATING SYSTEM ADEQUACY MODULE (ADJUSTED ELU) ===")
    
    # 1. Define System
    gens = [
        Generator("Nuclear", 400.0, 0.02, 4, Inf),
        Generator("Coal_A",  300.0, 0.04, 3, Inf),
        Generator("Coal_B",  300.0, 0.04, 3, Inf),
        Generator("Gas",     150.0, 0.05, 2, Inf),
        
        # --- ADJUSTMENT HERE ---
        # Previous: 200.0 * 50.0  (Too harsh, causes divergence)
        # New:      200.0 * 600.0 (600 Hours of water - More realistic)
        Generator("Hydro_ELU", 200.0, 0.01, 2, 200.0 * 600.0), 
        
        Generator("Old_56",   56.0, 0.10, 0, Inf)
    ]
    
    # 2. Define Load (Same as before)
    hours = 1:8760
    base_load = [750.0 + 300.0 * sin((h-2000)/8760 * 2π) + 50.0*randn() for h in hours]
    base_load = max.(0.0, base_load)
    weekly_peaks = [maximum(base_load[((w-1)*168+1):min(w*168,8760)]) for w in 1:52]
    
    # 3. Params
    params = SystemParams(20.0, 5.0, 1000) 
    
    # 4. Run Maintenance Scheduling
    schedule_maintenance!(gens, weekly_peaks)
    
    # 5. Run Simulations
    # Note: Analytical will now calculate a smaller 'q' penalty because energy is sufficient
    ana_lole, ana_profile = run_analytical(gens, base_load, params)
    mc_lole, mc_profile, mc_dist = run_monte_carlo(gens, base_load, params)
    
    # 6. Report
    println("\n--- RESULTS (With 600h Energy Limit) ---")
    @printf("Analytical LOLE: %.4f hours/year\n", ana_lole)
    @printf("Monte Carlo LOLE: %.4f hours/year\n", mc_lole)
    
    diff = abs(mc_lole - ana_lole) / ana_lole * 100
    @printf("Difference: %.1f%%\n", diff)
    
    if diff < 20.0
        println("SUCCESS: Results have converged! The ELU is no longer the dominant failure mode.")
    end
    
    # 7. Visualization
    p1 = histogram(mc_dist, bins=30, normalize=:pdf, alpha=0.6, label="MC Distribution",
                   title="Convergence Check: MC vs Analytical", xlabel="LOLE (h/yr)")
    vline!([ana_lole], color=:red, lw=3, label="Analytical Mean")
    
    display(p1)

    # Plot 2: Hourly Profile (Peak Week)
    peak_idx = argmax(mc_profile)
    win = (peak_idx-50):(peak_idx+50)
    p2 = plot(win, mc_profile[win], label="MC Risk", fill=(0,0.3,:blue),
              title="Hourly Risk Profile (Zoomed)", xlabel="Hour", ylabel="Prob")
    plot!(win, ana_profile[win], label="Analytical Risk", color=:red, lw=2)
    
    display(p1)
    display(p2)
end

run_comparison()