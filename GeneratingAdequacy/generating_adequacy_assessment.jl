using Plots
using Printf

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

struct Generator
    capacity::Float64
    for_rate::Float64 # Forced Outage Rate (q)
    name::String
end

struct COPT
    capacity_outage::Vector{Float64} # The X-axis (MW Outage)
    probability::Vector{Float64}     # The Y-axis (Probability)
end

# ==============================================================================
# 2. CORE ALGORITHM: DISCRETE CONVOLUTION (UNIT ADDITION)
# ==============================================================================

"""
    add_unit(copt::COPT, unit::Generator, step_size::Float64)

Adds a generator to the existing Capacity Outage Probability Table (COPT)
using the recursive convolution formula:
P_new(X) = P_old(X)*(1-q) + P_old(X-C)*q
"""
function add_unit(old_copt::COPT, unit::Generator, step_size::Float64)
    C = unit.capacity
    q = unit.for_rate
    p = 1.0 - q
    
    # 1. Determine the new maximum outage level
    # The new table must extend to (Old_Max + Unit_Capacity)
    max_outage_old = isempty(old_copt.capacity_outage) ? 0.0 : maximum(old_copt.capacity_outage)
    max_outage_new = max_outage_old + C
    
    # Create the new discrete states
    # We round the max size to ensure it fits the step grid
    num_steps = Int(ceil(max_outage_new / step_size)) + 1
    new_states = [i * step_size for i in 0:(num_steps-1)]
    new_probs = zeros(length(new_states))
    
    # 2. Helper to get probability from old table (handling boundary conditions)
    # Returns 0.0 if the requested outage state X doesn't exist in old table
    function get_old_prob(X_val)
        # Find index corresponding to X_val (approximate match for floats)
        idx = findfirst(x -> abs(x - X_val) < 1e-5, old_copt.capacity_outage)
        if idx === nothing
            return 0.0
        else
            return old_copt.probability[idx]
        end
    end

    # 3. Perform Convolution
    # P_new(X) = P_old(X) * p  +  P_old(X - C) * q
    
    # NOTE: Handling Rounding Logic (Slide 26-32)
    # If C is not a multiple of step_size, we must split the transition.
    # C_lower = floor(C), C_upper = ceil(C)
    # alpha = (C_upper - C) / step_size (Wait, slide says inverse ratio)
    
    # Let's implement the Exact Rounding Logic from the slides:
    # State 2 (C) is split into 2' (C_lower) and 2'' (C_upper)
    # Transition rate to C_lower is modified by alpha
    
    lower_idx = Int(floor(C / step_size))
    upper_idx = Int(ceil(C / step_size))
    C_lower = lower_idx * step_size
    C_upper = upper_idx * step_size
    
    if C_lower == C_upper
        # Case A: Exact Multiple (No Rounding needed)
        for (i, X) in enumerate(new_states)
            term1 = get_old_prob(X) * p             # Unit UP
            term2 = get_old_prob(X - C) * q         # Unit DOWN
            new_probs[i] = term1 + term2
        end
    else
        # Case B: Rounding Needed (Slide 30)
        # alpha is the weight for the UPPER state in terms of probability mass?
        # Slide 30: "Transition rates modified in inverse ratio of differences"
        # Let's use standard interpolation logic which matches the mean:
        # P(C_upper) * C_upper + P(C_lower) * C_lower = P_total * C
        
        # alpha = (C - C_lower) / (C_upper - C_lower)
        # This alpha goes to the UPPER state
        alpha = (C - C_lower) / step_size
        
        q_upper = q * alpha
        q_lower = q * (1.0 - alpha)
        
        println("  > Rounding Unit $(unit.name) ($(C) MW): Split into $(C_lower) MW and $(C_upper) MW")
        
        for (i, X) in enumerate(new_states)
            term1 = get_old_prob(X) * p                     # Unit UP
            term2 = get_old_prob(X - C_lower) * q_lower     # Unit DOWN (Rounded Low)
            term3 = get_old_prob(X - C_upper) * q_upper     # Unit DOWN (Rounded High)
            new_probs[i] = term1 + term2 + term3
        end
    end
    
    return COPT(new_states, new_probs)
end

# ==============================================================================
# 3. RELIABILITY INDICES CALCULATION
# ==============================================================================

function calculate_indices(copt::COPT, load_duration_curve::Vector{Float64})
    # LOLE: Loss of Load Expectation (Hours/Year)
    # EUE: Expected Unserved Energy (MWh/Year)
    
    lole = 0.0
    eue = 0.0
    
    # For every hour in the year (Load Duration Curve)
    for load in load_duration_curve
        # Probability of Failure = P(Available Capacity < Load)
        # Equivalently: P(Outage > Installed_Cap - Load)
        
        installed_cap = copt.capacity_outage[end] # Assuming max outage = total cap
        reserve = installed_cap - load
        
        # Sum probabilities where Outage > Reserve
        prob_loss_of_load = 0.0
        expected_shortage = 0.0
        
        for (i, outage) in enumerate(copt.capacity_outage)
            if outage > reserve
                prob = copt.probability[i]
                prob_loss_of_load += prob
                shortage_amount = outage - reserve
                expected_shortage += (shortage_amount * prob)
            end
        end
        
        lole += prob_loss_of_load
        eue += expected_shortage
    end
    
    return lole, eue
end

# ==============================================================================
# 4. MAIN EXECUTION
# ==============================================================================

# --- A. Define System ---
# Let's define a small test system
step_size = 10.0 # MW
gens = [
    Generator(50.0, 0.02, "Gen1"),
    Generator(50.0, 0.02, "Gen2"),
    Generator(56.0, 0.04, "Gen3_Odd"), # This will trigger rounding
    Generator(100.0, 0.05, "Gen4")
]

total_capacity = sum(g.capacity for g in gens)
println("Total Installed Capacity: $total_capacity MW")

# --- B. Build COPT ---
# Initialize with 0 MW outage having probability 1.0
current_copt = COPT([0.0], [1.0])

println("\n--- Building COPT ---")
for g in gens
    global current_copt
    current_copt = add_unit(current_copt, g, step_size)
    println("Added $(g.name): Table size is now $(length(current_copt.capacity_outage)) states.")
end

# --- C. Define Load ---
# Simple Load Duration Curve: 8760 hours
# Peak 200 MW, Min 100 MW, Linear decrease
hours = 1:8760
peak_load = 200.0
min_load = 100.0
slope = (peak_load - min_load) / 8760
load_curve = [peak_load - slope * (h-1) for h in hours]

# --- D. Calculate Indices ---
println("\n--- Calculating Indices ---")
lole, eue = calculate_indices(current_copt, load_curve)

@printf("LOLE (Loss of Load Expectation): %.4f hours/year\n", lole)
@printf("EUE (Expected Unserved Energy):  %.4f MWh/year\n", eue)

# --- E. Visualization ---
# 1. Plot COPT (Cumulative Probability)
# We usually plot P(Outage >= X) vs X
cumulative_probs = zeros(length(current_copt.probability))
cum_sum = 0.0
# Sum from right to left
for i in length(current_copt.probability):-1:1
    global cum_sum
    cum_sum += current_copt.probability[i]
    cumulative_probs[i] = cum_sum
end

p1 = plot(current_copt.capacity_outage, log10.(cumulative_probs), 
    title="Capacity Outage Probability Table (COPT)",
    xlabel="Capacity Outage (MW)", ylabel="Log10(Cumulative Probability)",
    label="P(Outage >= X)", linewidth=2, color=:blue, legend=:topright)

# 2. Plot Load vs Capacity
p2 = plot(hours, load_curve, fill=(0, 0.2, :green), label="Load Duration Curve",
     xlabel="Hours", ylabel="MW", title="System Adequacy Overview")
hline!([total_capacity], label="Total Installed Capacity", color=:red, linestyle=:dash)

display(p1)
display(p2)
