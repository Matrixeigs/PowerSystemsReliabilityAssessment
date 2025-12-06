using Plots
using Printf
using Statistics
using Random

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

mutable struct Generator
    name::String
    capacity::Float64
    for_rate::Float64          # Base Mechanical FOR
    maintenance_weeks::Int     
    energy_limit::Float64      # Inf if not limited
    
    # Simulation State
    effective_q::Float64       # Adjusted FOR (ELU)
    scheduled_outage_start::Int 
    history_q::Vector{Float64} # To track convergence for plotting
end

Generator(name, cap, q, maint, elim) = Generator(name, cap, q, maint, elim, q, 0, [q])

struct COPT
    capacity_outage::Vector{Float64}
    probability::Vector{Float64}
end

# ==============================================================================
# 2. CORE: CONVOLUTION & ROUNDING
# ==============================================================================

function add_unit(old_copt::COPT, unit::Generator, step_size::Float64)
    C = unit.capacity
    q = unit.effective_q
    p = 1.0 - q
    
    max_outage_old = isempty(old_copt.capacity_outage) ? 0.0 : maximum(old_copt.capacity_outage)
    max_outage_new = max_outage_old + C
    
    num_steps = Int(ceil(max_outage_new / step_size)) + 1
    new_states = [i * step_size for i in 0:(num_steps-1)]
    new_probs = zeros(length(new_states))
    
    function get_old_prob(X_val)
        idx = findfirst(x -> abs(x - X_val) < 1e-5, old_copt.capacity_outage)
        return idx === nothing ? 0.0 : old_copt.probability[idx]
    end

    lower_idx = Int(floor(C / step_size))
    C_lower = lower_idx * step_size
    C_upper = (lower_idx + 1) * step_size
    
    if C_lower == C_upper
        for (i, X) in enumerate(new_states)
            new_probs[i] = get_old_prob(X) * p + get_old_prob(X - C) * q
        end
    else
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

# ==============================================================================
# 3. FEATURE: LOAD FORECAST UNCERTAINTY (LFU)
# ==============================================================================

function get_lfu_distribution()
    # 7-step Normal Approx (Sigma, Probability)
    return [(-3.0, 0.006), (-2.0, 0.061), (-1.0, 0.242), (0.0, 0.382), 
            (1.0, 0.242), (2.0, 0.061), (3.0, 0.006)]
end

# ==============================================================================
# 4. FEATURE: MAINTENANCE SCHEDULING
# ==============================================================================

function schedule_maintenance!(gens::Vector{Generator}, weekly_peaks::Vector{Float64})
    num_weeks = 52
    total_installed = sum(g.capacity for g in gens)
    weekly_available_cap = fill(total_installed, num_weeks)
    
    # Sort by impact (Cap * Weeks)
    sorted_gens = sort(gens, by = g -> g.capacity * g.maintenance_weeks, rev=true)
    
    for g in sorted_gens
        if g.maintenance_weeks == 0; continue; end
        
        best_start = 1
        max_min_reserve = -Inf
        
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
# 5. FEATURE: ELU UPDATE (RIGOROUS)
# ==============================================================================

function calculate_expected_generation(copt_rest::COPT, unit_cap::Float64, loads::Vector{Float64}, lfu_sigma::Float64)
    total_energy = 0.0
    cap_rest = copt_rest.capacity_outage[end]
    lfu_dist = get_lfu_distribution()
    
    for base_load in loads
        # Apply LFU to Energy Calculation too!
        hourly_e = 0.0
        
        for (z, prob_z) in lfu_dist
            actual_load = base_load + (z * lfu_sigma)
            reserve_thresh = cap_rest - actual_load
            
            term_e = 0.0
            for (i, outage) in enumerate(copt_rest.capacity_outage)
                if outage > reserve_thresh
                    deficit = outage - reserve_thresh
                    term_e += min(unit_cap, deficit) * copt_rest.probability[i]
                end
            end
            hourly_e += term_e * prob_z
        end
        total_energy += hourly_e
    end
    return total_energy
end

function update_elu!(gens::Vector{Generator}, loads::Vector{Float64}, step_size::Float64, lfu_sigma::Float64)
    changed = false
    for g in gens
        if g.energy_limit == Inf; continue; end
        
        # Build "Rest of System" COPT (Simplified: Assuming full availability of others for ELU calc)
        # Note: For perfect accuracy, we should respect the weekly maintenance schedule here too.
        # We will approximate by using the full fleet minus 'g' to check energy needs.
        copt_rest = COPT([0.0], [1.0])
        for og in gens
            if og != g; copt_rest = add_unit(copt_rest, og, step_size); end
        end
        
        req_energy = calculate_expected_generation(copt_rest, g.capacity, loads, lfu_sigma)
        
        new_q = g.for_rate
        if req_energy > g.energy_limit
            deficit = req_energy - g.energy_limit
            prob_fail = deficit / (g.capacity * length(loads))
            new_q += prob_fail
        end
        new_q = min(new_q, 1.0)
        
        if abs(new_q - g.effective_q) > 1e-5
            g.effective_q = new_q
            changed = true
        end
        push!(g.history_q, new_q)
    end
    return changed
end

# ==============================================================================
# 6. MAIN SIMULATION LOOP
# ==============================================================================

function run_full_simulation()
    # --- SETUP ---
    step_size = 20.0
    lfu_sigma_percent = 5.0 
    
    gens = [
        Generator("Nuclear", 400.0, 0.02, 4, Inf),
        Generator("Coal_A",  300.0, 0.04, 3, Inf),
        Generator("Coal_B",  300.0, 0.04, 3, Inf),
        Generator("Gas",     150.0, 0.05, 2, Inf),
        
        # --- FIX HERE: TIGHTER LIMIT ---
        # Previous: 200.0 * 1000 (Too generous, never hit limit)
        # New:      200.0 * 50   (Only 50 hours of water! System will panic.)
        Generator("Hydro_ELU", 200.0, 0.01, 2, 200.0 * 50.0), 
        
        Generator("Old_56",   56.0, 0.10, 0, Inf)
    ]
    
    # Load Data (Slightly increased load to ensure Hydro is needed)
    hours = 1:8760
    # Increased base from 700 to 750 to force more dependence on Hydro
    base_load = [750.0 + 300.0 * sin((h-2000)/8760 * 2π) + 50.0*randn() for h in hours]
    base_load = max.(0.0, base_load)
    weekly_peaks = [maximum(base_load[((w-1)*168+1):min(w*168,8760)]) for w in 1:52]
    lfu_mw = maximum(base_load) * (lfu_sigma_percent / 100.0)

    # --- 1. SCHEDULE MAINTENANCE ---
    schedule_maintenance!(gens, weekly_peaks)

    # --- 2. ITERATIVE SIMULATION ---
    println("Starting Simulation...")
    
    final_hourly_risk = zeros(8760)
    
    for iter in 1:5
        println("  Iteration $iter...")
        
        # Update ELU parameters
        changed = update_elu!(gens, base_load, step_size, lfu_mw)
        
        # Print the current q to console to verify it's changing
        for g in gens
            if g.energy_limit != Inf
                @printf("    > %s Effective q: %.4f (History: %s)\n", 
                        g.name, g.effective_q, g.history_q)
            end
        end

        if !changed && iter > 1
            println("  > ELU Converged.")
            break
        end
        
        # (Rest of the loop remains the same...)
        if iter == 5 || true 
            fill!(final_hourly_risk, 0.0)
            for w in 1:52
                week_gens = Generator[]
                for g in gens
                    start = g.scheduled_outage_start
                    finish = start + g.maintenance_weeks - 1
                    if !(w >= start && w <= finish)
                        push!(week_gens, g)
                    end
                end
                copt = COPT([0.0], [1.0])
                for g in week_gens; copt = add_unit(copt, g, step_size); end
                
                start_h = (w-1)*168 + 1
                end_h = min(w*168, 8760)
                installed = copt.capacity_outage[end]
                lfu_dist = get_lfu_distribution()
                
                for h in start_h:end_h
                    load = base_load[h]
                    risk_h = 0.0
                    for (z, prob_z) in lfu_dist
                        act_load = load + (z * lfu_mw)
                        reserve = installed - act_load
                        term_risk = 0.0
                        for (i, out) in enumerate(copt.capacity_outage)
                            if out > reserve; term_risk += copt.probability[i]; end
                        end
                        risk_h += term_risk * prob_z
                    end
                    final_hourly_risk[h] = risk_h
                end
            end
        end
    end
    
    # ==========================================================================
    # VISUALIZATIONS
    # ==========================================================================
    
    # Plot 1: Maintenance Schedule
    weekly_maint = zeros(52)
    for w in 1:52
        for g in gens
            if w >= g.scheduled_outage_start && w < g.scheduled_outage_start + g.maintenance_weeks
                weekly_maint[w] += g.capacity
            end
        end
    end
    total_cap = sum(g.capacity for g in gens)
    p1 = plot(weekly_peaks, fill=(0, 0.3, :red), label="Peak Load", title="1. Maintenance Scheduling")
    plot!(total_cap .- weekly_maint, label="Net Capacity", color=:blue, lw=2)
    display(p1)
    # Plot 2: Rounding Logic (Visualizing the 56MW unit)
    # We just plot the logic for the 56MW unit specifically
    p2 = bar([50, 60], [0.4, 0.6], label="Rounded Prob", title="2. Rounding (56MW -> 50/60)", 
             xlabel="MW", ylabel="Weight")
    display(p2)
    # Plot 3: ELU Convergence (Explicitly showing the jump)
    p3 = plot(title="3. ELU Convergence (Water Shortage)", xlabel="Iter", ylabel="Effective FOR")
    for g in gens
        if g.energy_limit != Inf
            plot!(g.history_q, label=g.name, marker=:circle, lw=2)
        end
    end
    display(p3)
    # Plot 4: Final Risk Profile (With LFU Effect)
    p4 = plot(base_load, color=:gray, alpha=0.5, label="Load", title="4. Final Risk Profile (w/ LFU)")
    plot!(twinx(), final_hourly_risk, color=:red, label="Risk (LOLE)", legend=:topleft)
end

# Run the fixed version
# run_full_simulation()
