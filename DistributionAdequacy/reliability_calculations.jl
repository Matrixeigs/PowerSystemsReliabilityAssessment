# =============================================================================
# Distribution System and Station Adequacy Assessment - Numerical Examples
# Based on the provided PDF content
# =============================================================================

using Printf

# =============================================================================
# PART 1: DISTRIBUTION SYSTEM ADEQUACY (Cases 1, 2, 3)
# =============================================================================

println("================================================================")
println("PART 1: DISTRIBUTION SYSTEM ADEQUACY ASSESSMENT")
println("================================================================")

# --- Data Definitions ---
# Component Data
# Main Feeder sections and Lateral sections are defined by failure rates and repair times
# Note: The PDF implies specific lengths or aggregate failure rates for sections.
# Based on Slide 19 calculations:
# Load Point A (Main + Lat A): lambda contribution seems to be 1.35? 
# Let's reverse engineer the exact input parameters from the solved example in Slide 19.
# Slide 19 says:
#   A: lambda=1.35, r=1.55 -> U = 2.0925 (approx 2.10)
#   B: lambda=1.10, r=2.05 -> U = 2.255 (approx 2.25)
#   C: lambda=0.85, r=2.05 -> U = 1.7425 (approx 1.75)

# Let's define the Load Points directly with these calculated aggregate parameters for Case 1
struct LoadPoint
    name::String
    customers::Int
    lambda::Float64 # failures/yr
    r::Float64      # average outage duration (hr)
    U::Float64      # annual outage time (hr/yr)
end

# Customer counts
N_A = 250
N_B = 100
N_C = 50
N_Total = 400

# --- Case 1: Fused Laterals, No Alternate Supply ---
println("\n--- Case 1: Fused Laterals, No Alternate Supply ---")

# Data from Slide 19
lp_A_c1 = LoadPoint("A", N_A, 1.35, 1.55, 1.35 * 1.55)
lp_B_c1 = LoadPoint("B", N_B, 1.10, 2.05, 1.10 * 2.05)
lp_C_c1 = LoadPoint("C", N_C, 0.85, 2.05, 0.85 * 2.05)

load_points_c1 = [lp_A_c1, lp_B_c1, lp_C_c1]

function calculate_system_indices(lps::Vector{LoadPoint}, total_cust::Int)
    sum_lambda_N = sum(lp.lambda * lp.customers for lp in lps)
    sum_U_N = sum(lp.U * lp.customers for lp in lps)
    
    saifi = sum_lambda_N / total_cust
    saidi = sum_U_N / total_cust
    caidi = saidi / saifi
    asai = (total_cust * 8760 - sum_U_N) / (total_cust * 8760)
    
    return saifi, saidi, caidi, asai
end

saifi_1, saidi_1, caidi_1, asai_1 = calculate_system_indices(load_points_c1, N_Total)

@printf("SAIFI: %.2f interruptions/system customer/yr\n", saifi_1)
@printf("SAIDI: %.2f hr/system customer/yr\n", saidi_1)
@printf("CAIDI: %.2f hr/customer interruption\n", caidi_1)
@printf("ASAI : %.6f\n", asai_1)


# --- Case 2: Alternate Supply Available ---
println("\n--- Case 2: Alternate Supply Available (Switching time = 1.0 hr) ---")
println("Note: Failure rates (lambda) remain the same. Repair times (r) decrease for Main Feeder faults due to switching.")

# In Case 2, faults on the main feeder can be restored by switching (1 hr) instead of repair (3 hr).
# This reduces the average outage duration 'r', and thus 'U'.
# We don't have the raw component breakdown to recalculate exactly, but we can infer the improvement.
# However, for the purpose of this script, we will use the logic described:
# "Load point failure rates are not affected... SAIDI decreases."

# Let's assume hypothetical improved 'r' values based on switching logic
# If Main Feeder fault (which affects everyone) is now 1hr instead of 3hr.
# Let's assume Main Feeder lambda was 0.1 * Length.
# Without exact component lengths, we'll print the qualitative result or placeholders.
println("Calculating with assumed improvement in restoration times...")

# Hypothetical improvement: Assume 50% of outage time was due to Main Feeder repair (3hr), now becomes 1hr.
# This is an estimation for demonstration as exact component lengths weren't in the text summary.
lp_A_c2 = LoadPoint("A", N_A, 1.35, 1.05, 1.35 * 1.05) # Improved r
lp_B_c2 = LoadPoint("B", N_B, 1.10, 1.55, 1.10 * 1.55) # Improved r
lp_C_c2 = LoadPoint("C", N_C, 0.85, 1.55, 0.85 * 1.55) # Improved r

load_points_c2 = [lp_A_c2, lp_B_c2, lp_C_c2]
saifi_2, saidi_2, caidi_2, asai_2 = calculate_system_indices(load_points_c2, N_Total)

@printf("SAIFI: %.2f (Unchanged from Case 1)\n", saifi_2)
@printf("SAIDI: %.2f (Improved)\n", saidi_2)


# --- Case 3: Solidly Connected (No Fuses) ---
println("\n--- Case 3: Solidly Connected (No Fuses) ---")
println("In this case, any fault on any lateral trips the main breaker.")
println("All load points experience the SUM of all failure rates.")

# Total system lambda = sum of all component lambdas
# From Case 1, we can estimate individual contributions roughly, but let's use the logic:
# Lambda_System = Lambda_Main + Lambda_LatA + Lambda_LatB + Lambda_LatC
# Let's assume the total failure rate is the sum of the max seen in Case 1 (approx 1.6 or higher).

# Using values typical for this comparison:
lambda_total = 1.60 # Estimated sum
r_avg = 2.0         # Average repair time

lp_A_c3 = LoadPoint("A", N_A, lambda_total, r_avg, lambda_total * r_avg)
lp_B_c3 = LoadPoint("B", N_B, lambda_total, r_avg, lambda_total * r_avg)
lp_C_c3 = LoadPoint("C", N_C, lambda_total, r_avg, lambda_total * r_avg)

load_points_c3 = [lp_A_c3, lp_B_c3, lp_C_c3]
saifi_3, saidi_3, caidi_3, asai_3 = calculate_system_indices(load_points_c3, N_Total)

@printf("SAIFI: %.2f (Significantly Higher)\n", saifi_3)
@printf("SAIDI: %.2f (Significantly Higher)\n", saidi_3)


# =============================================================================
# PART 2: STATION RELIABILITY (Simple Configuration)
# =============================================================================

println("\n================================================================")
println("PART 2: STATION RELIABILITY ASSESSMENT (Case a)")
println("================================================================")

# Configuration (a): Source -> Breaker -> Transformer -> Bus -> Load
# This is a Series System.

# Component Data (Table 10)
# Lambda (f/yr), Repair Time (hr)
# Note: Breaker has Active and Passive. For series calc, we sum them roughly for continuity.

lambda_source = 0.0   # Assumed reliable
lambda_breaker_act = 0.006
lambda_breaker_pas = 0.004
lambda_breaker = lambda_breaker_act + lambda_breaker_pas
r_breaker = 10.0

lambda_trans = 0.010
r_trans = 100.0

lambda_bus = 0.001
r_bus = 2.0

# Calculation for Load Point
# For a series system:
# Lambda_Load = Sum(Lambda_Components)
# U_Load = Sum(Lambda_Component * r_Component)
# r_Load = U_Load / Lambda_Load

lambda_station_load = lambda_breaker + lambda_trans + lambda_bus
U_station_load = (lambda_breaker * r_breaker) + (lambda_trans * r_trans) + (lambda_bus * r_bus)
r_station_load = U_station_load / lambda_station_load

println("Configuration: Source -> Breaker -> Transformer -> Bus -> Load")
println("Calculation Method: Series System (Approximate)")

@printf("Load Point Failure Rate (lambda): %.4f f/yr\n", lambda_station_load)
@printf("Load Point Outage Duration (r)  : %.2f hr\n", r_station_load)
@printf("Load Point Unavailability (U)   : %.4f hr/yr\n", U_station_load)

println("\nDone.")