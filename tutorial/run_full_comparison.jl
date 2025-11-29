include("PowerSystemAdequacy.jl")
using .PowerSystemAdequacy

# 1. Define System
#    (ID, Cap, MTTF, MTTR)
gen_data = [
    (1, 400.0, 1100.0, 50.0),
    (2, 400.0, 1100.0, 50.0),
    (3, 300.0, 1200.0, 60.0),
    (4, 300.0, 1200.0, 60.0),
    (5, 150.0, 900.0,  40.0),
    (6, 150.0, 900.0,  40.0),
    (7, 50.0,  500.0,  20.0), # Small unit
    (8, 50.0,  500.0,  20.0)
]
gens = [Generator(id, cap, mttf, mttr) for (id, cap, mttf, mttr) in gen_data]

# 2. Load Model
hours = 1:8760
# Peak 1600 MW
raw_load = 1100.0 .+ 500.0 .* sin.((hours .- 2000) ./ 8760 .* 2π) .+ 100.0 .* randn(8760)
load_model = LoadModel(max.(0.0, raw_load))

# 3. Run All Methods
#    A. Analytical (Exact)
res_ana = run_analytical(gens, load_model, step_size=10.0)

#    B. Non-Sequential (Fast MC) - 5000 iterations
res_nsmc = run_non_sequential_mc(gens, load_model, 5000)

#    C. Sequential (Detailed MC) - 500 years
res_smc = run_sequential_mc(gens, load_model, 500)

# 4. Compare
compare_results([res_ana, res_nsmc, res_smc])
