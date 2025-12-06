using Printf
using DataFrames

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

"""
Generator with Frequency & Duration Parameters
"""
struct GeneratorFD
    name::String
    capacity::Float64   # MW
    mtbf::Float64       # Mean Time Between Failures (Hours) = 1/lambda
    mttr::Float64       # Mean Time To Repair (Hours) = 1/mu
    
    # Derived Parameters (Calculated in constructor)
    lambda::Float64     # Failures/Year
    mu::Float64         # Repairs/Year
    p::Float64          # Availability
    q::Float64          # Unavailability (FOR)
    
    function GeneratorFD(name, cap, mtbf_h, mttr_h)
        # Convert hours to rates (per year)
        # Note: 8760 hours/year is standard
        lam = 8760.0 / mtbf_h
        mu_rate = 8760.0 / mttr_h
        
        # Steady State Probabilities
        q = lam / (lam + mu_rate)
        p = 1.0 - q
        
        new(name, cap, mtbf_h, mttr_h, lam, mu_rate, p, q)
    end
end

"""
Capacity Outage Probability & Frequency Table (COPT)
"""
struct COPT
    outage_levels::Vector{Float64} # MW
    cum_prob::Vector{Float64}      # P(Outage >= X)
    cum_freq::Vector{Float64}      # F(Outage >= X) [occ/yr]
end

# ==============================================================================
# 2. RECURSIVE ALGORITHM (With Educational Reporting)
# ==============================================================================

"""
Adds a unit to the COPT and prints the step-by-step calculation logic.
"""
function add_unit_educational!(copt::COPT, unit::GeneratorFD)
    println("\n" * "="^80)
    @printf("ADDING UNIT: %s (Cap: %.0f MW, q: %.4f, λ: %.2f f/yr)\n", 
            unit.name, unit.capacity, unit.q, unit.lambda)
    println("="^80)
    
    # 1. Prepare New Arrays
    # The new max outage is OldMax + UnitCapacity
    current_max = maximum(copt.outage_levels)
    new_max = current_max + unit.capacity
    
    # Create a dense grid of outage levels (assuming integer steps for simplicity)
    # In real tools, we use sparse steps, but for education, 1MW or 10MW steps are fine.
    # Here we define the significant states based on existing + new combinations.
    
    # For simplicity in this demo, we assume integer MW steps
    step_size = 1.0 
    new_levels = collect(0.0:step_size:new_max)
    
    new_probs = zeros(length(new_levels))
    new_freqs = zeros(length(new_levels))
    
    # Helper to get Old P and F (Interpolated or Nearest)
    function get_old_metrics(x)
        if x < 0
            # Boundary Condition: Outage < 0 is impossible to be "less than", 
            # so P(Outage >= negative) = 1.0 (System always has outage >= -5)
            # F(Outage >= negative) = 0.0 (Never crosses from < -5 to >= -5)
            return 1.0, 0.0
        end
        
        # Find index
        idx = findfirst(v -> v >= x, copt.outage_levels)
        if isnothing(idx)
            return 0.0, 0.0 # Outage > Max Possible
        end
        
        # Exact match check (simplified)
        if abs(copt.outage_levels[idx] - x) < 1e-5
            return copt.cum_prob[idx], copt.cum_freq[idx]
        else
            # If x is between levels, P(>=x) is the same as P(>=next_level) 
            # strictly speaking for discrete, but let's assume step function
            return copt.cum_prob[idx], copt.cum_freq[idx]
        end
    end

    # 2. Reporting Header
    println("\n--- Step-by-Step Calculation (Selected States) ---")
    @printf("%-10s | %-25s | %-35s\n", "Outage(X)", "Probability P(X)", "Frequency F(X)")
    println("-"^80)
    
    # 3. Perform Recursion
    C = unit.capacity
    p, q = unit.p, unit.q
    lam = unit.lambda
    
    for (i, X) in enumerate(new_levels)
        # Get Old Values
        P_old_X, F_old_X = get_old_metrics(X)
        P_old_X_C, F_old_X_C = get_old_metrics(X - C)
        
        # --- A. Probability Recursion ---
        # P_new(X) = p * P_old(X) + q * P_old(X - C)
        P_new = (p * P_old_X) + (q * P_old_X_C)
        
        # --- B. Frequency Recursion ---
        # F_new(X) = p * F_old(X) + q * F_old(X - C) + λ * p * [P_old(X - C) - P_old(X)]
        term1 = p * F_old_X
        term2 = q * F_old_X_C
        term3 = lam * p * (P_old_X_C - P_old_X)
        
        F_new = term1 + term2 + term3
        
        new_probs[i] = P_new
        new_freqs[i] = F_new
        
        # Educational Print for significant states
        # We print if X is a multiple of Capacity or near 0
        if mod(X, max(1, C)) == 0 || X == 0
            # Format the math string
            p_math = @sprintf("%.2f*%.4f + %.2f*%.4f", p, P_old_X, q, P_old_X_C)
            f_math = @sprintf("T1:%.2f + T2:%.2f + T3:%.2f", term1, term2, term3)
            
            @printf("%-10.0f | %-25s = %-6.4f | %-35s = %-6.4f\n", 
                    X, "Calc..", P_new, f_math, F_new)
        end
    end
    
    println("-"^80)
    println("Note on Frequency Terms:")
    println("  T1 = p * F_old(X)        (System fluctuates while Unit is UP)")
    println("  T2 = q * F_old(X-C)      (System fluctuates while Unit is DOWN)")
    println("  T3 = λ * p * [P(X-C)-P(X)] (Unit fails, pushing system into outage)")
    
    return COPT(new_levels, new_probs, new_freqs)
end

# ==============================================================================
# 3. EVALUATION FUNCTIONS
# ==============================================================================

function evaluate_risk(copt::COPT, peak_load::Float64, installed_cap::Float64)
    # Reserve = Installed - Peak
    # System Fails if Outage > Reserve
    reserve = installed_cap - peak_load
    
    # Find index in COPT where Outage > Reserve
    # Since COPT is "Cumulative >= X", we need the first X that is > Reserve.
    # Actually, Loss of Load occurs if Available < Load <=> (Cap - Outage) < Load
    # <=> Outage > (Cap - Load) <=> Outage > Reserve.
    
    # We need P(Outage > Reserve). 
    # In discrete steps, P(Outage > R) is P(Outage >= Next_Step_Above_R)
    
    idx = findfirst(x -> x > reserve, copt.outage_levels)
    
    if isnothing(idx)
        return 0.0, 0.0, 0.0
    end
    
    LOLE_prob = copt.cum_prob[idx]
    LOLF_occ  = copt.cum_freq[idx]
    
    # Calculate Duration
    # Duration = Probability / Frequency (Hours/Year / Occ/Year -> Hours/Occ)
    # Note: LOLE is usually in Hours/Year. Here P is probability [0-1].
    # So LOLE_hours = P * 8760
    
    LOLE_hours = LOLE_prob * 8760.0
    LOLD_hours = (LOLF_occ > 0) ? (LOLE_hours / LOLF_occ) : 0.0
    
    return LOLE_hours, LOLF_occ, LOLD_hours
end

# ==============================================================================
# 4. MAIN EXECUTION (CASE STUDY)
# ==============================================================================

function run_educational_demo()
    println("### GENERATING SYSTEM ADEQUACY: FREQUENCY & DURATION DEMO ###")
    
    # 1. Initialize Empty System (0 MW Outage, P=1.0, F=0.0)
    # We define a base COPT with 0 capacity
    copt = COPT([0.0], [1.0], [0.0])
    
    total_cap = 0.0
    
    # 2. Define Units (Matches your PPT Example)
    # Unit 1: 16MW, lambda=2, mu=98
    # MTTF = 8760/2 = 4380h, MTTR = 8760/98 = 89.38h
    u1 = GeneratorFD("Unit 1", 16.0, 4380.0, 89.39) 
    
    # Unit 2: Identical
    u2 = GeneratorFD("Unit 2", 16.0, 4380.0, 89.39)
    
    # 3. Add Unit 1
    copt = add_unit_educational!(copt, u1)
    total_cap += u1.capacity
    
    # 4. Add Unit 2
    copt = add_unit_educational!(copt, u2)
    total_cap += u2.capacity
    
    # 5. Risk Assessment
    peak_load = 20.0 # MW
    println("\n" * "="^80)
    println("RISK ASSESSMENT (Peak Load = $peak_load MW, Installed = $total_cap MW)")
    println("Reserve = $(total_cap - peak_load) MW")
    println("="^80)
    
    lole, lolf, lold = evaluate_risk(copt, peak_load, total_cap)
    
    @printf("LOLE (Expectation):  %8.4f hours/year\n", lole)
    @printf("LOLF (Frequency):    %8.4f occ/year\n", lolf)
    @printf("LOLD (Mean Duration):%8.4f hours/occurrence\n", lold)
    println("="^80)
end

# Run the demo
run_educational_demo()
