# =============================================================================
# Reliability Assessment for Distribution Systems Using Linear Programming
# Non-Simulation-Based Approach
# =============================================================================
#
# Reference:
# G. Muñoz-Delgado, J. Contreras, and J. M. Arroyo,
# "Reliability Assessment for Distribution Optimization Models:
# A Non-Simulation-Based Linear Programming Approach,"
# IEEE Transactions on Smart Grid, vol. 9, no. 4, pp. 3048-3059, 2018.
#
# Converted from MATLAB to Julia
# =============================================================================

using JuMP
using HiGHS
using SparseArrays
using Printf
using LinearAlgebra
using Plots
using Statistics
# =============================================================================
# DATA STRUCTURES
# =============================================================================

"""
    DistributionSystem

Data structure for distribution system network.

Fields:
- bus: Bus data [BUS_I, NC, PD, BUS_E]
- branch: Branch data [F_BUS, T_BUS, LENGTH, TAU_RP, TAU_SW, F_BUS_E, T_BUS_E]
- gen: Generator/substation data [GEN_BUS, GEN_BUS_E]
- lambda: Failure rate per unit length (failures/year/km)
"""
struct DistributionSystem
    bus::Matrix{Float64}
    branch::Matrix{Float64}
    gen::Matrix{Float64}
    lambda::Float64
end

"""
    ReliabilityIndices

Reliability indices for distribution system.
"""
struct ReliabilityIndices
    SAIFI::Vector{Float64}  # System Average Interruption Frequency Index per customer
    SAIDI::Vector{Float64}  # System Average Interruption Duration Index per customer
    ASAI::Vector{Float64}   # Average Service Availability Index per customer
    EENS::Float64           # Expected Energy Not Supplied (MWh/year)
    CIFs::Vector{Float64}   # Customer Interruption Frequency
    CIDs::Vector{Float64}   # Customer Interruption Duration
end

# Bus column indices
const BUS_I = 1
const NC = 2
const PD = 3
const BUS_E = 4

# Branch column indices
const F_BUS = 1
const T_BUS = 2
const LENGTH = 3
const TAU_RP = 4
const TAU_SW = 5
const F_BUS_E = 6
const T_BUS_E = 7

# Generator column indices
const GEN_BUS = 1
const GEN_BUS_E = 2

# =============================================================================
# TEST CASE DATA
# =============================================================================

"""
    create_test_case_33()

Create IEEE 33-bus radial distribution test system for reliability assessment.
"""
function create_test_case_33()
    # Bus data: [BUS_I, NC (customers), PD (kW), BUS_E]
    bus = Float64[
        1   0     0      1;
        2   100   100    2;
        3   90    90     3;
        4   120   120    4;
        5   60    60     5;
        6   60    60     6;
        7   200   200    7;
        8   200   200    8;
        9   60    60     9;
        10  60    60     10;
        11  45    45     11;
        12  60    60     12;
        13  60    60     13;
        14  120   120    14;
        15  60    60     15;
        16  60    60     16;
        17  60    90     17;
        18  90    90     18;
        19  90    90     19;
        20  90    90     20;
        21  90    90     21;
        22  90    90     22;
        23  90    90     23;
        24  420   420    24;
        25  420   420    25;
        26  60    60     26;
        27  60    60     27;
        28  60    120    28;
        29  120   200    29;
        30  200   150    30;
        31  150   210    31;
        32  210   60     32;
        33  60    60     33
    ]

    # Branch data: [F_BUS, T_BUS, LENGTH(km), TAU_RP(hrs), TAU_SW(hrs), F_BUS_E, T_BUS_E]
    branch = Float64[
        1   2   0.5   3.0   1.0   1   2;
        2   3   0.5   3.0   1.0   2   3;
        3   4   0.5   3.0   1.0   3   4;
        4   5   0.5   3.0   1.0   4   5;
        5   6   0.5   3.0   1.0   5   6;
        6   7   0.5   3.0   1.0   6   7;
        7   8   1.0   3.0   1.0   7   8;
        8   9   0.5   3.0   1.0   8   9;
        9   10  0.5   3.0   1.0   9   10;
        10  11  0.5   3.0   1.0   10  11;
        11  12  0.5   3.0   1.0   11  12;
        12  13  0.5   3.0   1.0   12  13;
        13  14  1.0   3.0   1.0   13  14;
        14  15  0.5   3.0   1.0   14  15;
        15  16  0.5   3.0   1.0   15  16;
        16  17  0.5   3.0   1.0   16  17;
        17  18  0.5   3.0   1.0   17  18;
        2   19  0.5   3.0   1.0   2   19;
        19  20  0.5   3.0   1.0   19  20;
        20  21  0.5   3.0   1.0   20  21;
        21  22  0.5   3.0   1.0   21  22;
        3   23  0.5   3.0   1.0   3   23;
        23  24  0.5   3.0   1.0   23  24;
        24  25  0.5   3.0   1.0   24  25;
        6   26  0.5   3.0   1.0   6   26;
        26  27  0.5   3.0   1.0   26  27;
        27  28  1.0   3.0   1.0   27  28;
        28  29  0.5   3.0   1.0   28  29;
        29  30  0.5   3.0   1.0   29  30;
        30  31  0.5   3.0   1.0   30  31;
        31  32  0.5   3.0   1.0   31  32;
        32  33  1.0   3.0   1.0   32  33
    ]

    # Generator/substation data: [GEN_BUS, GEN_BUS_E]
    gen = Float64[
        1  1
    ]

    # Failure rate (failures/year/km)
    lambda = 0.2

    return DistributionSystem(bus, branch, gen, lambda)
end

# =============================================================================
# MAIN RELIABILITY ASSESSMENT FUNCTION
# =============================================================================

"""
    reliability_assessment_linear_programming(mpc::DistributionSystem)

Perform reliability assessment using non-simulation-based linear programming approach.

Returns:
- ReliabilityIndices: Calculated reliability metrics
- computation_time: Elapsed time for computation
"""
function reliability_assessment_linear_programming(mpc::DistributionSystem)
    println("\n" * "="^70)
    println("RELIABILITY ASSESSMENT - LINEAR PROGRAMMING APPROACH")
    println("="^70)

    time_start = time()

    bus = copy(mpc.bus)
    branch = copy(mpc.branch)
    gen = copy(mpc.gen)
    lambda = mpc.lambda

    # Data sizes
    nb = size(bus, 1)
    nl = size(branch, 1)
    ng = size(gen, 1)
    nd = count(bus[:, PD] .> 0)

    @printf("System size: %d buses, %d branches, %d generators, %d load buses\n", nb, nl, ng, nd)

    # Validation
    if nl + ng != nb
        error("The size of input information is incorrect!")
    end

    # Convert external to internal bus numbering
    i2e = Int.(bus[:, BUS_I])
    e2i = zeros(Int, maximum(i2e))
    for (i, ext) in enumerate(i2e)
        e2i[ext] = i
    end

    bus[:, BUS_I] = collect(1:nb)

    branch[:, F_BUS] = [e2i[Int(b)] for b in branch[:, F_BUS_E]]
    branch[:, T_BUS] = [e2i[Int(b)] for b in branch[:, T_BUS_E]]

    gen[:, GEN_BUS] = [e2i[Int(b)] for b in gen[:, GEN_BUS_E]]
    gen_bus = Int.(gen[:, GEN_BUS])

    # Detect feeders and substations
    index_ss = Int[]
    index_child = Int[]

    for g in 1:ng
        index_f = findall(branch[:, F_BUS] .== gen_bus[g])
        index_child_g = Int.(branch[index_f, T_BUS])

        index_t = findall(branch[:, T_BUS] .== gen_bus[g])
        append!(index_child, index_child_g)
        append!(index_child, Int.(branch[index_t, F_BUS]))
        append!(index_ss, index_f)
        append!(index_ss, index_t)
    end

    nf = length(index_ss)

    @printf("Number of feeders: %d\n", nf)

    # Construct connection matrices
    Cg = sparse(gen_bus, 1:ng, ones(ng), nb, ng)

    load_buses = findall(bus[:, PD] .> 0)
    Cd = sparse(load_buses, 1:nd, ones(nd), nb, nd)
    load_bus = Int.(bus[load_buses, BUS_I])

    # Connection matrix for network
    f = Int.(branch[:, F_BUS])
    t = Int.(branch[:, T_BUS])
    i = vcat(1:nl, 1:nl)
    Cft = sparse(i, vcat(f, t), vcat(ones(nl), -ones(nl)), nl, nb)

    # Branch parameters
    lambda_ij = branch[:, LENGTH] * lambda
    tau_RP_ij = branch[:, TAU_RP]
    tau_SW_ij = branch[:, TAU_SW]

    # Define optimization problem dimensions
    FIJ = 0
    FJI = FIJ + nl
    FG = FJI + nl
    nx = FG + ng

    # Setup linear programming model for path generation
    println("\nGenerating paths for each load bus...")

    fd = zeros(nb)
    Aeq_span = spzeros(nb, nx)
    beq_span = -fd

    Aeq_span[:, FIJ+1:FJI] = Cft'
    Aeq_span[:, FJI+1:FG] = -Cft'
    Aeq_span[:, FG+1:nx] = -Cg

    cobj = ones(nx)
    cobj[FG+1:nx] .= 0

    # Generate paths for each load bus
    path_set = zeros(nd, nl)

    for d in 1:nd
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, 0 <= x[1:nx] <= 1)
        @objective(model, Min, cobj' * x)

        # Equality constraints
        beq_span_d = copy(beq_span)
        beq_span_d[load_bus[d]] = -1.0

        @constraint(model, Aeq_span * x .== beq_span_d)

        optimize!(model)

        if termination_status(model) == MOI.OPTIMAL
            x_sol = value.(x)
            path_set[d, :] = x_sol[FIJ+1:FJI] + x_sol[FJI+1:FG]
        else
            @warn "No optimal solution found for load bus $d"
        end
    end

    println("Path generation complete.")

    # Analyze reliability indices
    println("\nCalculating reliability indices...")

    # (1) Expected rate of nodal repair-and-switching interruptions
    Ns_RP = path_set * lambda_ij

    # (2) Expected rate of switching-only interruptions
    Nij_CB = zeros(nl)
    path_set_T = copy(path_set)

    for f in 1:nf
        path_set_T[:, index_ss[f]] .= 0
        if index_child[f] <= nd
            path_set_T[index_child[f], index_ss[f]] = 1
        end
    end

    Nij_CB = path_set_T * lambda_ij
    Ns_T = zeros(nf)

    for f in 1:nf
        Ns_T[f] = path_set[:, index_ss[f]]' * Nij_CB
    end

    Ns_SW = zeros(nd)
    for d in 1:nd
        feeder_indices = findall(path_set[d, index_ss] .> 0)
        if !isempty(feeder_indices)
            Ns_SW[d] = Ns_T[feeder_indices[1]] - Ns_RP[d]
        end
    end

    # (3) Expected duration of nodal repair-and-switching interruptions
    Ds_RP = path_set * (lambda_ij .* tau_RP_ij)

    # (4) Expected duration of switching-only interruptions
    Ds_SW_bar = path_set_T * (lambda_ij .* tau_SW_ij)
    Ds_T = zeros(nf)

    for f in 1:nf
        Ds_T[f] = path_set[:, index_ss[f]]' * Ds_SW_bar
    end

    Ds_WW = path_set * (lambda_ij .* tau_SW_ij)
    Ds_SW = zeros(nd)

    for d in 1:nd
        feeder_indices = findall(path_set[d, index_ss] .> 0)
        if !isempty(feeder_indices)
            Ds_SW[d] = Ds_T[feeder_indices[1]] - Ds_WW[d]
        end
    end

    # Derive reliability indices
    CIFs = Ns_RP + Ns_SW  # Customer Interruption Frequency
    CIDs = Ds_RP + Ds_SW  # Customer Interruption Duration

    consumer = bus[load_bus, NC]
    total_consumers = sum(consumer)

    SAIFI = (consumer .* CIFs) / total_consumers
    SAIDI = (consumer .* CIDs) / total_consumers
    ASAI = 1.0 .- SAIDI / 8760.0

    # Expected Energy Not Supplied (EENS)
    mu_b = [0.70, 0.83, 1.00]  # Load factors
    Delta_b = [2000.0, 5760.0, 1000.0]  # Time periods
    EENS = sum((CIDs .* bus[load_bus, PD]) * mu_b' * Delta_b) / 8760.0 / 1000.0

    computation_time = time() - time_start

    # Print summary results
    println("\n" * "="^70)
    println("RELIABILITY ASSESSMENT RESULTS")
    println("="^70)
    @printf("System SAIFI: %.4f interruptions/customer/year\n", sum(SAIFI))
    @printf("System SAIDI: %.4f hours/customer/year\n", sum(SAIDI))
    @printf("System ASAI:  %.6f\n", mean(ASAI))
    @printf("System EENS:  %.4f MWh/year\n", EENS)
    @printf("\nComputation time: %.3f seconds\n", computation_time)
    println("="^70)

    indices = ReliabilityIndices(SAIFI, SAIDI, ASAI, EENS, CIFs, CIDs)

    return indices, computation_time
end

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

"""
    plot_reliability_indices(indices::ReliabilityIndices, output_dir::String="figures")

Generate plots for reliability indices.
"""
function plot_reliability_indices(indices::ReliabilityIndices, output_dir::String="figures")
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    nd = length(indices.CIFs)

    # Plot 1: Customer Interruption Frequency
    p1 = bar(1:nd, indices.CIFs, xlabel="Load Bus", ylabel="CIF (interruptions/year)",
             title="Customer Interruption Frequency", color=:steelblue, legend=false)
    savefig(p1, joinpath(output_dir, "fig_lp_cif.png"))

    # Plot 2: Customer Interruption Duration
    p2 = bar(1:nd, indices.CIDs, xlabel="Load Bus", ylabel="CID (hours/year)",
             title="Customer Interruption Duration", color=:coral, legend=false)
    savefig(p2, joinpath(output_dir, "fig_lp_cid.png"))

    # Plot 3: SAIFI per customer
    p3 = bar(1:nd, indices.SAIFI, xlabel="Load Bus", ylabel="SAIFI (int/cust/yr)",
             title="System Average Interruption Frequency Index", color=:mediumseagreen,
             legend=false)
    savefig(p3, joinpath(output_dir, "fig_lp_saifi.png"))

    # Plot 4: SAIDI per customer
    p4 = bar(1:nd, indices.SAIDI, xlabel="Load Bus", ylabel="SAIDI (hr/cust/yr)",
             title="System Average Interruption Duration Index", color=:orchid,
             legend=false)
    savefig(p4, joinpath(output_dir, "fig_lp_saidi.png"))

    println("\nGenerated 4 reliability index plots in $output_dir/")
end

# =============================================================================
# MAIN EXECUTION
# =============================================================================

"""
    main()

Main execution function.
"""
function main()
    println("\n" * "="^80)
    println(" LINEAR PROGRAMMING APPROACH FOR DISTRIBUTION RELIABILITY")
    println("="^80)

    # Create test case
    println("\nLoading test case...")
    mpc = create_test_case_33()

    # Run reliability assessment
    indices, comp_time = reliability_assessment_linear_programming(mpc)

    # Generate visualizations
    println("\nGenerating visualizations...")
    plot_reliability_indices(indices, "figures")

    println("\n" * "="^80)
    println(" ANALYSIS COMPLETE")
    println("="^80)
    println()

    return indices
end

# Run the main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
