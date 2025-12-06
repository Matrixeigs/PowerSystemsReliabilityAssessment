using Plots
using LinearAlgebra
using Random
using Statistics

# Set a fixed seed for reproducibility
Random.seed!(42)

# ==============================================================================
# PART 1 & 2: DEFINITIONS AND RATES
# ==============================================================================
println("--- PART 1 & 2: MTTF, MTTR, and RATES ---")

# 1. Define Parameters (Based on typical values from your MATLAB scripts)
# Let's take a Generator as an example
MTTF_val = 1000.0  # Mean Time To Failure (Hours)
MTTR_val = 50.0    # Mean Time To Repair (Hours)

# 2. Calculate Rates
# The fundamental link: Rate = 1 / Mean Time
λ = 1 / MTTF_val  # Failure Rate (failures per hour)
μ = 1 / MTTR_val  # Repair Rate (repairs per hour)

println("Given Parameters:")
println("  MTTF: $MTTF_val hours")
println("  MTTR: $MTTR_val hours")
println("Derived Rates:")
println("  Failure Rate (λ): $(round(λ, digits=5)) /hour")
println("  Repair Rate (μ):  $(round(μ, digits=5)) /hour")
println("--------------------------------------------------\n")

# ==============================================================================
# PART 3: WHY CONSTANT RATE -> EXPONENTIAL CDF
# ==============================================================================
println("--- PART 3: SIMULATION PROOF OF EXPONENTIAL DISTRIBUTION ---")
println("Simulating 10,000 components with a CONSTANT probability of failure...")

# Simulation parameters
N_samples = 10000
dt = 1.0 # Time step (1 hour)
max_time = 5000

# We simulate the "Coin Flip" logic:
# Every hour, check if it fails with probability (λ * dt).
# This assumes NO memory.
failure_times = Float64[]

for i in 1:N_samples
    t = 0.0
    while true
        # The "Constant Rate" Assumption:
        # Random check against constant threshold
        if rand() < (λ * dt) 
            push!(failure_times, t)
            break
        end
        t += dt
        if t > max_time break end
    end
end

# Visualization
histogram(failure_times, normalize=:pdf, label="Simulated (Constant Rate)", 
          title="Why Constant Rate = Exponential Distribution",
          xlabel="Time to Failure (Hours)", ylabel="Probability Density",
          alpha=0.6, color=:blue, bins=50)

# Overlay the Theoretical Exponential PDF: f(t) = λ * exp(-λt)
t_theory = 0:10:max_time
pdf_theory = λ .* exp.(-λ .* t_theory)
p3 = plot!(t_theory, pdf_theory, label="Theoretical Exp PDF", 
           linewidth=3, color=:red, linestyle=:dash)

display(p3)
println(">> Plot 1 Generated: Histogram matches the theoretical curve.")
println("--------------------------------------------------\n")

# ==============================================================================
# PART 4: SINGLE COMPONENT MARKOV PROCESS (Analytical vs Monte Carlo)
# ==============================================================================
println("--- PART 4: SINGLE COMPONENT MARKOV CHAIN ---")

# 1. Define the Transition Matrix P for dt = 1 hour
# State 0: Working (Up), State 1: Failed (Down)
# P = [ p00  p01 ]
#     [ p10  p11 ]

# Using the exact exponential formula derived in the chat
p01 = 1 - exp(-λ * dt) # Prob Up -> Down
p10 = 1 - exp(-μ * dt) # Prob Down -> Up
p00 = 1 - p01          # Stay Up
p11 = 1 - p10          # Stay Down

P = [p00 p01; 
     p10 p11]

println("Transition Matrix P (1-step):")
display(P)

# 2. Analytical Evolution (Matrix Multiplication)
steps = 200
prob_down_analytical = zeros(steps)
current_state_prob = [1.0 0.0] # Start 100% Working [Prob_Up, Prob_Down]

for t in 1:steps
    global current_state_prob
    # Markov Equation: π(t+1) = π(t) * P
    current_state_prob = current_state_prob * P 
    prob_down_analytical[t] = current_state_prob[2]
end

# 3. Monte Carlo Simulation (Sequential Sampling like MATLAB)
# We simulate one component moving through time
mc_state_history = zeros(steps)
current_state = 0 # 0 = Up, 1 = Down

for t in 1:steps
    global current_state
    r = rand()
    if current_state == 0 # Currently Up
        if r < p01
            current_state = 1 # Fail
        end
    else # Currently Down
        if r < p10
            current_state = 0 # Repair
        end
    end
    mc_state_history[t] = current_state
end

# Visualization
p4 = plot(1:steps, prob_down_analytical, label="Analytical Prob(Down)", 
     title="Markov Process: Single Component", xlabel="Time (Hours)", ylabel="Probability / State",
     linewidth=3, color=:blue)
plot!(1:steps, mc_state_history, label="Monte Carlo Realization (0/1)", 
      seriestype=:steppost, color=:orange, alpha=0.5, linestyle=:dot)
# Add Steady State Line
steady_state_unavailability = MTTR_val / (MTTF_val + MTTR_val)
hline!([steady_state_unavailability], label="Steady State Limit", color=:green)

display(p4)
println(">> Plot 2 Generated: Convergence of probability vs Random state jumps.")
println("--------------------------------------------------\n")

# ==============================================================================
# PART 5: MULTI-COMPONENT SYSTEM (Sequential Monte Carlo)
# ==============================================================================
println("--- PART 5: MULTI-COMPONENT SYSTEM SIMULATION ---")
println("Simulating a system with 5 Generators over 1 year (8760 hours)...")

# System Definition
num_gens = 5
# Let's give them slightly different reliability data
gen_mttfs = [1000.0, 1200.0, 800.0, 1500.0, 2000.0]
gen_mttrs = [50.0,   60.0,   40.0,  20.0,   100.0]
gen_capacities = [100, 100, 50, 200, 150] # MW

# Pre-calculate transition probabilities for each generator
# Matrix size: num_gens x 2 (Col 1: P_fail, Col 2: P_repair)
trans_probs = zeros(num_gens, 2)
for i in 1:num_gens
    trans_probs[i, 1] = 1 - exp(-(1/gen_mttfs[i])) # p01
    trans_probs[i, 2] = 1 - exp(-(1/gen_mttrs[i])) # p10
end

# Simulation Loop (Sequential Monte Carlo)
hours_year = 1000 # Shortened for visualization clarity (usually 8760)
system_available_capacity = zeros(hours_year)
gen_states = zeros(Int, num_gens) # All start at 0 (Up)

for t in 1:hours_year
    # Update every generator
    for i in 1:num_gens
        r = rand()
        if gen_states[i] == 0 # If Up
            if r < trans_probs[i, 1]
                gen_states[i] = 1 # Fail
            end
        else # If Down
            if r < trans_probs[i, 2]
                gen_states[i] = 0 # Repair
            end
        end
    end
    
    # Calculate System Capacity for this hour
    current_cap = 0.0
    for i in 1:num_gens
        if gen_states[i] == 0
            current_cap += gen_capacities[i]
        end
    end
    system_available_capacity[t] = current_cap
end

# Visualization
total_cap = sum(gen_capacities)
p5 = plot(1:hours_year, system_available_capacity, 
     title="System Capacity Availability (Markov Process)",
     xlabel="Time (Hours)", ylabel="Available MW",
     label="Available Capacity", fill=(0, 0.2, :blue), size=(800, 400))
hline!([total_cap], label="Max Capacity", color=:red, linestyle=:dash)

display(p5)
println(">> Plot 3 Generated: System Capacity fluctuating over time.")
println("   Notice how individual Markov processes combine to form system state.")
