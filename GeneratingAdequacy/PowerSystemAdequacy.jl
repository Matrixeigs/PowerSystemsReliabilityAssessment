module PowerSystemAdequacy

using Statistics
using Random
using Printf
using Plots

export Generator, LoadModel, ReliabilityResult, 
       run_analytical, run_non_sequential_mc, run_sequential_mc, 
       compare_results

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

"""
    Generator
defined by Capacity and Reliability Parameters (MTTF/MTTR).
"""
struct Generator
    id::Int
    capacity::Float64   # MW
    mttf::Float64       # Hours
    mttr::Float64       # Hours
    
    # Derived Parameters
    lambda::Float64     # Failure Rate (1/MTTF)
    mu::Float64         # Repair Rate (1/MTTR)
    for_rate::Float64   # Forced Outage Rate (Unavailability)
end

function Generator(id::Int, capacity::Float64, mttf::Float64, mttr::Float64)
    λ = 1.0 / mttf
    μ = 1.0 / mttr
    u = λ / (λ + μ)
    return Generator(id, capacity, mttf, mttr, λ, μ, u)
end

struct LoadModel
    hourly_load::Vector{Float64}
    peak_load::Float64
    function LoadModel(hourly_load::Vector{Float64})
        new(hourly_load, maximum(hourly_load))
    end
end

struct ReliabilityResult
    method::String
    lole_hours_yr::Float64
    eue_mwh_yr::Float64
    computation_time::Float64
    convergence_history::Vector{Float64} # Only for MC
end

# Helper for Analytical Method
struct COPT
    capacity_outage::Vector{Float64}
    probability::Vector{Float64}
end

# ==============================================================================
# 2. ENGINE A: ANALYTICAL (Convolution)
#    - Exact calculation (within rounding error)
#    - No sampling noise
# ==============================================================================

function add_unit_convolution(old_copt::COPT, unit::Generator, step_size::Float64)
    C = unit.capacity
    q = unit.for_rate
    p = 1.0 - q
    
    # Determine new table size
    max_outage_old = isempty(old_copt.capacity_outage) ? 0.0 : maximum(old_copt.capacity_outage)
    max_outage_new = max_outage_old + C
    num_steps = Int(ceil(max_outage_new / step_size)) + 1
    
    new_states = [i * step_size for i in 0:(num_steps-1)]
    new_probs = zeros(length(new_states))
    
    # Helper to find probability in old table
    function get_prob(X_val)
        # Simple index lookup assuming sorted, fixed steps
        idx = Int(round(X_val / step_size)) + 1
        if idx < 1 || idx > length(old_copt.probability)
            return 0.0
        end
        return old_copt.probability[idx]
    end

    # Rounding Logic (Splitting probability between steps if Capacity isn't a multiple)
    lower_idx = Int(floor(C / step_size))
    C_lower = lower_idx * step_size
    C_upper = (lower_idx + 1) * step_size
    
    if abs(C - C_lower) < 1e-5 # Exact Match
        for (i, X) in enumerate(new_states)
            new_probs[i] = get_prob(X) * p + get_prob(X - C) * q
        end
    else # Interpolation (Rounding)
        alpha = (C - C_lower) / step_size
        q_upper = q * alpha
        q_lower = q * (1.0 - alpha)
        for (i, X) in enumerate(new_states)
            new_probs[i] = get_prob(X) * p + 
                           get_prob(X - C_lower) * q_lower + 
                           get_prob(X - C_upper) * q_upper
        end
    end
    
    return COPT(new_states, new_probs)
end

function run_analytical(gens::Vector{Generator}, load::LoadModel; step_size::Float64=10.0)
    t_start = time()
    println("--- Running Analytical Method (Convolution) ---")
    
    # 1. Build System COPT
    copt = COPT([0.0], [1.0])
    for g in gens
        copt = add_unit_convolution(copt, g, step_size)
    end
    
    # 2. Evaluate Risk against Load Profile
    total_installed = sum(g.capacity for g in gens)
    lole = 0.0
    eue = 0.0
    
    # Pre-calculate Cumulative Probability for speed
    # P(Outage >= X)
    cum_prob = reverse(cumsum(reverse(copt.probability)))
    
    for load_mw in load.hourly_load
        reserve = total_installed - load_mw
        
        # We need P(Outage > Reserve)
        # Find index in COPT where Outage > Reserve
        # Since states are 0, step, 2*step...
        # We need index > Reserve / step
        
        idx = Int(floor(reserve / step_size)) + 2 # +1 for 0-based, +1 for > logic
        
        if idx <= length(cum_prob) && idx >= 1
            prob_loss = cum_prob[idx]
            lole += prob_loss
            
            # EUE requires iterating the tail
            # Simplified EUE: sum((Outage - Reserve) * prob) for relevant states
            for k in idx:length(copt.capacity_outage)
                outage = copt.capacity_outage[k]
                eue += (outage - reserve) * copt.probability[k]
            end
        elseif idx < 1
            # Reserve is negative (Load > Installed), 100% loss
            lole += 1.0
            eue += (load_mw - total_installed) # Base deficit
            # Plus average outage
            avg_outage = sum(copt.capacity_outage .* copt.probability)
            eue += avg_outage
        end
    end
    
    return ReliabilityResult("Analytical", lole, eue, time() - t_start, Float64[])
end

# ==============================================================================
# 3. ENGINE B: NON-SEQUENTIAL MC
# ==============================================================================

function run_non_sequential_mc(gens::Vector{Generator}, load::LoadModel, iterations::Int)
    t_start = time()
    println("--- Running Non-Sequential MC ($iterations iters) ---")
    
    cum_lole = 0.0
    cum_eue = 0.0
    history = Float64[]
    
    total_hours = length(load.hourly_load)
    
    for i in 1:iterations
        # Sample State
        cap_avail = 0.0
        for g in gens
            if rand() >= g.for_rate
                cap_avail += g.capacity
            end
        end
        
        # Evaluate Year
        iter_lole = 0.0
        iter_eue = 0.0
        for h in 1:total_hours
            if cap_avail < load.hourly_load[h]
                deficit = load.hourly_load[h] - cap_avail
                iter_lole += 1.0
                iter_eue += deficit
            end
        end
        
        cum_lole += iter_lole
        cum_eue += iter_eue
        
        if i % 100 == 0
            push!(history, cum_lole / i)
        end
    end
    
    return ReliabilityResult("Non-Sequential MC", cum_lole/iterations, cum_eue/iterations, time()-t_start, history)
end

# ==============================================================================
# 4. ENGINE C: SEQUENTIAL MC
# ==============================================================================

function run_sequential_mc(gens::Vector{Generator}, load::LoadModel, years::Int)
    t_start = time()
    println("--- Running Sequential MC ($years years) ---")
    
    n_gens = length(gens)
    total_hours = length(load.hourly_load)
    
    # State: Time until next transition
    # status: true=UP, false=DOWN
    status = trues(n_gens)
    ttf = [-log(rand())/g.lambda for g in gens]
    
    cum_lole = 0.0
    cum_eue = 0.0
    history = Float64[]
    
    for y in 1:years
        year_lole = 0.0
        year_eue = 0.0
        
        for h in 1:total_hours
            # Update States
            cap_avail = 0.0
            for i in 1:n_gens
                g = gens[i]
                ttf[i] -= 1.0
                while ttf[i] <= 0
                    if status[i] # Was UP, fails
                        status[i] = false
                        ttf[i] += -log(rand())/g.mu
                    else # Was DOWN, repairs
                        status[i] = true
                        ttf[i] += -log(rand())/g.lambda
                    end
                end
                if status[i]; cap_avail += g.capacity; end
            end
            
            # Evaluate
            if cap_avail < load.hourly_load[h]
                deficit = load.hourly_load[h] - cap_avail
                year_lole += 1.0
                year_eue += deficit
            end
        end
        
        cum_lole += year_lole
        cum_eue += year_eue
        
        if y % 10 == 0
            push!(history, cum_lole / y)
        end
    end
    
    return ReliabilityResult("Sequential MC", cum_lole/years, cum_eue/years, time()-t_start, history)
end

# ==============================================================================
# 5. COMPARISON TOOLS
# ==============================================================================

function compare_results(results::Vector{ReliabilityResult})
    println("\n==========================================")
    println("       METHOD COMPARISON SUMMARY")
    println("==========================================")
    @printf "%-20s | %-10s | %-10s | %-10s\n" "Method" "LOLE(h/yr)" "EUE(MWh)" "Time(s)"
    println("-"^60)
    
    for r in results
        @printf "%-20s | %-10.4f | %-10.2f | %-10.4f\n" r.method r.lole_hours_yr r.eue_mwh_yr r.computation_time
    end
    println("-"^60)
    
    # Plot Convergence for MC methods
    p = plot(title="Monte Carlo Convergence", xlabel="Samples (x10/x100)", ylabel="LOLE", lw=2)
    for r in results
        if !isempty(r.convergence_history)
            plot!(r.convergence_history, label=r.method)
        else
            # Plot Analytical as a straight line
            hline!([r.lole_hours_yr], label="Analytical (Benchmark)", linestyle=:dash, color=:black)
        end
    end
    display(p)
end

end # module
