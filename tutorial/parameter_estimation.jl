using Plots
using Statistics
using Random

Random.seed!(123) # Fix seed for reproducible "field data"

# ==============================================================================
# PART 1: VISUALIZING FIELD DATA COLLECTION (THE LOG)
# ==============================================================================
println("--- GENERATING FIELD DATA LOGS ---")

# 1. Simulate a "Real" Machine History
# We assume hidden true physics, but we will "measure" them from the plot
true_mttf = 200.0 
true_mttr = 50.0
n_events = 6

# Generate alternating Up (1) and Down (0) intervals
# We use Exponential distribution for the simulation source
up_durations = -true_mttf .* log.(rand(n_events))
down_durations = -true_mttr .* log.(rand(n_events))

# Construct the Timeline for Plotting
time_points = [0.0]
state_values = [1.0] # Start Up
current_t = 0.0

# We will store data to annotate the plot
annotation_points = [] 

for i in 1:n_events
    # 1. Run period (Up)
    global current_t += up_durations[i]
    push!(time_points, current_t, current_t)
    push!(state_values, 1.0, 0.0) # Drop to 0
    push!(annotation_points, (current_t - up_durations[i]/2, 1.05, "TTF_$i"))
    
    # 2. Repair period (Down)
    global current_t += down_durations[i]
    push!(time_points, current_t, current_t)
    push!(state_values, 0.0, 1.0) # Rise to 1
    push!(annotation_points, (current_t - down_durations[i]/2, 0.1, "TTR_$i"))
end

# --- FIGURE 1: THE OPERATIONAL TIMELINE ---
p1 = plot(time_points, state_values, 
    seriestype=:steppre, 
    title="Figure 1: Field Data Collection (Operational Log)",
    label="Machine State (1=Up, 0=Down)",
    xlabel="Time (Hours)", ylabel="State",
    linewidth=2, color=:green,
    fill=(0, 0.2, :green), # Fill area to show "Up" time clearly
    xlims=(0, time_points[end]), ylims=(0, 1.3),
    size=(800, 400),
    legend=:topright
)

# Add Red shading for Down times (manually for visual clarity)
# And add annotations
for i in 1:n_events
    # Annotate TTF (Green zones)
    annotate!(annotation_points[i*2-1]...)
    # Annotate TTR (Red zones)
    annotate!(annotation_points[i*2]...)
end

display(p1)
println(">> Figure 1 Generated: This represents the raw logs engineers analyze.")
println("   - Green Areas = Time To Failure (TTF)")
println("   - Gaps = Time To Repair (TTR)")
println("--------------------------------------------------\n")


# ==============================================================================
# PART 2: CALCULATING RATES FROM DATA (CONVERGENCE)
# ==============================================================================
println("--- CALCULATING PARAMETERS FROM DATA ---")

# Now we simulate a long history to show how we derive the Rate
# As we collect more data samples, our calculated Rate approaches the truth.

n_long_run = 1000
long_up_times = -true_mttf .* log.(rand(n_long_run))

# Arrays to store the "Running Calculation"
calculated_mttf = zeros(n_long_run)
calculated_lambda = zeros(n_long_run)

current_sum = 0.0
for i in 1:n_long_run
    global current_sum += long_up_times[i]
    
    # FORMULA 1: Mean Time = Total Time / Number of Events
    avg_ttf = current_sum / i
    calculated_mttf[i] = avg_ttf
    
    # FORMULA 2: Failure Rate = 1 / Mean Time
    calculated_lambda[i] = 1 / avg_ttf
end

# --- FIGURE 2: ESTIMATION CONVERGENCE ---
p2 = plot(1:n_long_run, calculated_lambda, 
    title="Figure 2: Deriving Failure Rate from Field Data",
    xlabel="Number of Observed Failures (Sample Size)", 
    ylabel="Calculated Failure Rate (λ)",
    label="Estimated λ = (N / Σ TTF)",
    linewidth=2, color=:blue,
    size=(800, 400)
)

# Add the "True" value line
hline!([1/true_mttf], label="True Theoretical Rate", color=:red, linestyle=:dash, linewidth=2)

display(p2)

# Final Calculations Display
est_mttf = mean(up_durations)
est_mttr = mean(down_durations)
est_lambda = 1 / est_mttf
est_mu = 1 / est_mttr

println("CALCULATION RESULTS (Based on Figure 1 Data):")
println("1. Sum of Up Times (Green)   = $(round(sum(up_durations), digits=2)) hours")
println("2. Number of Failures        = $n_events")
println("3. Calculated MTTF           = $(round(est_mttf, digits=2)) hours")
println("4. DERIVED FAILURE RATE (λ)  = 1 / MTTF = $(round(est_lambda, digits=5)) failures/hr")
println("   (True value was $(1/true_mttf))")
println("")
println("5. Sum of Down Times (Gaps)  = $(round(sum(down_durations), digits=2)) hours")
println("6. Number of Repairs         = $n_events")
println("7. Calculated MTTR           = $(round(est_mttr, digits=2)) hours")
println("8. DERIVED REPAIR RATE (μ)   = 1 / MTTR = $(round(est_mu, digits=5)) repairs/hr")
