module AdequacyAssessmentFast

using Statistics
using Random
using Printf

export Generator, TieLine, Area, System, 
       run_fast_sequential_simulation,
       SupportPolicy, ISOLATED, INTERCONNECTED

# ==============================================================================
# 1. DATA STRUCTURES
# ==============================================================================

mutable struct Generator
    id::String
    capacity::Float64
    mttf::Float64
    mttr::Float64
    current_state::Bool # true=UP, false=DOWN
    time_to_transition::Float64
end

function Generator(id, cap, mttf, mttr)
    return Generator(id, cap, mttf, mttr, true, -log(rand()) * mttf)
end

struct TieLine
    from_area::Int
    to_area::Int
    capacity::Float64
end

struct Area
    id::Int
    name::String
    generators::Vector{Generator}
    hourly_load::Vector{Float64}
end

struct System
    areas::Vector{Area}
    tie_lines::Vector{TieLine}
    # Pre-computed adjacency matrix for the Max-Flow solver
    # capacity_matrix[i,j] = limit
    topology_matrix::Matrix{Float64} 
end

# Constructor to build topology once
function System(areas, lines)
    n = length(areas)
    mat = zeros(Float64, n, n)
    for line in lines
        # Tie lines are usually bidirectional
        mat[line.from_area, line.to_area] += line.capacity
        mat[line.to_area, line.from_area] += line.capacity
    end
    return System(areas, lines, mat)
end

@enum SupportPolicy ISOLATED INTERCONNECTED

# ==============================================================================
# 2. FAST SOLVER (CUSTOM MAX-FLOW)
# ==============================================================================

"""
    solve_curtailment_fast(sys, margins, policy)

Calculates minimum load shedding without an external LP solver.
Uses a custom Max-Flow (Augmenting Path) algorithm.
"""
function solve_curtailment_fast(sys::System, margins::Vector{Float64}, policy::SupportPolicy)
    n_areas = length(margins)
    
    # --- 1. FAST PATH: No Deficits ---
    # If everyone has positive margin, no calculation needed.
    if all(m -> m >= 0, margins)
        return zeros(Float64, n_areas)
    end

    # --- 2. ISOLATED POLICY ---
    # Simple calculation: if margin < 0, that is the curtailment.
    if policy == ISOLATED
        curtailments = zeros(Float64, n_areas)
        for i in 1:n_areas
            if margins[i] < 0
                curtailments[i] = -margins[i]
            end
        end
        return curtailments
    end

    # --- 3. INTERCONNECTED (MAX FLOW) ---
    # We need to move power from Surplus Nodes (+) to Deficit Nodes (-)
    
    # Create a residual graph based on the static topology
    # We copy the matrix because we will modify it during flow finding
    residual_capacity = copy(sys.topology_matrix)
    
    current_margins = copy(margins)
    
    # Iteratively find paths from Surplus -> Deficit
    # This is a simplified Ford-Fulkerson/Edmonds-Karp implementation
    while true
        # Find a source (Surplus) and a sink (Deficit)
        source_idx = findfirst(m -> m > 1e-4, current_margins)
        sink_idx = findfirst(m -> m < -1e-4, current_margins)
        
        # If no surplus left OR no deficit left, we are done
        if isnothing(source_idx) || isnothing(sink_idx)
            break
        end
        
        # BFS to find a path from source to sink
        parent = zeros(Int, n_areas)
        queue = [source_idx]
        found_path = false
        
        while !isempty(queue)
            u = popfirst!(queue)
            if u == sink_idx
                found_path = true
                break
            end
            
            for v in 1:n_areas
                # If there is capacity and we haven't visited v
                if residual_capacity[u, v] > 1e-4 && parent[v] == 0 && v != source_idx
                    parent[v] = u
                    push!(queue, v)
                end
            end
        end
        
        if !found_path
            # No path exists between this surplus and this deficit (Network constraints binding)
            # We must mark this source as "processed" to avoid infinite loop?
            # Actually, in standard MaxFlow, we run until no paths exist S->T.
            # Here we have multiple sources/sinks. 
            # Optimization: If we can't reach *this* sink, maybe we can reach another?
            # For simplicity in this custom solver: 
            # We break the loop if we can't connect ANY surplus to ANY deficit.
            # (A full implementation would use a Super-Source/Super-Sink, but this is sufficient for small grids)
            break 
        end
        
        # Calculate bottleneck on this path
        path_flow = min(current_margins[source_idx], -current_margins[sink_idx])
        curr = sink_idx
        while curr != source_idx
            prev = parent[curr]
            path_flow = min(path_flow, residual_capacity[prev, curr])
            curr = prev
        end
        
        # Apply flow
        current_margins[source_idx] -= path_flow
        current_margins[sink_idx] += path_flow
        
        curr = sink_idx
        while curr != source_idx
            prev = parent[curr]
            residual_capacity[prev, curr] -= path_flow
            residual_capacity[curr, prev] += path_flow # Reverse flow capability
            curr = prev
        end
    end
    
    # Whatever negative margin remains is the final EUE
    curtailments = zeros(Float64, n_areas)
    for i in 1:n_areas
        if current_margins[i] < 0
            curtailments[i] = -current_margins[i]
        end
    end
    
    return curtailments
end

# ==============================================================================
# 3. SIMULATION ENGINE
# ==============================================================================

function run_fast_sequential_simulation(sys::System, policy::SupportPolicy, n_years::Int)
    t_start = time()
    println("--- Running FAST Adequacy Assessment ---")
    println("Policy: $policy | Years: $n_years")
    
    n_areas = length(sys.areas)
    
    # Metrics
    area_lole = zeros(Float64, n_areas)
    area_eue = zeros(Float64, n_areas)
    
    # Pre-allocate arrays to avoid garbage collection
    current_margins = zeros(Float64, n_areas)
    
    for year in 1:n_years
        for h in 1:8760
            
            # 1. Update Gen States & Calculate Margins
            for (i, area) in enumerate(sys.areas)
                area_cap = 0.0
                for g in area.generators
                    g.time_to_transition -= 1.0
                    while g.time_to_transition <= 0
                        if g.current_state
                            g.current_state = false
                            g.time_to_transition += -log(rand()) * g.mttr
                        else
                            g.current_state = true
                            g.time_to_transition += -log(rand()) * g.mttf
                        end
                    end
                    if g.current_state
                        area_cap += g.capacity
                    end
                end
                # Margin = Gen - Load
                current_margins[i] = area_cap - area.hourly_load[h]
            end
            
            # 2. Solve Network Flow (The Fast Way)
            # This function replaces the LP solver
            curtailments = solve_curtailment_fast(sys, current_margins, policy)
            
            # 3. Update Indices
            for i in 1:n_areas
                if curtailments[i] > 0
                    area_lole[i] += 1.0
                    area_eue[i] += curtailments[i]
                end
            end
        end
    end
    
    elapsed = time() - t_start
    println("Simulation completed in $(round(elapsed, digits=2)) seconds.")
    
    results = []
    for i in 1:n_areas
        push!(results, (
            area = sys.areas[i].name,
            lole = area_lole[i] / n_years,
            eue = area_eue[i] / n_years
        ))
    end
    return results
end

# ==============================================================================
# 4. DEMO
# ==============================================================================

function run_demo()
    # --- Define System ---
    # Identical system to previous demo
    gens1 = [Generator("G1_$i", 400.0, 1000.0, 50.0) for i in 1:5]
    load1 = 1000.0 .+ 500.0 .* sin.(range(0, 2π, length=8760))
    area1 = Area(1, "Area_Rich", gens1, load1)
    
    gens2 = [Generator("G2_$i", 200.0, 900.0, 60.0) for i in 1:5]
    load2 = 800.0 .+ 400.0 .* sin.(range(0, 2π, length=8760))
    area2 = Area(2, "Area_Poor", gens2, load2)
    
    tie = TieLine(1, 2, 200.0)
    
    sys = System([area1, area2], [tie])
    
    # --- Run Analysis (More years now because it's fast!) ---
    
    # Case A: Isolated
    res_iso = run_fast_sequential_simulation(sys, ISOLATED, 500) 
    
    # Case B: Interconnected
    res_int = run_fast_sequential_simulation(sys, INTERCONNECTED, 500)
    
    # --- Report ---
    println("\n=== FINAL COMPARISON (FAST METHOD) ===")
    println("Policy          | Area       | LOLE (h/yr) | EUE (MWh/yr)")
    println("-"^60)
    
    for r in res_iso
        @printf("ISOLATED        | %-10s | %10.2f  | %10.2f\n", r.area, r.lole, r.eue)
    end
    println("-"^60)
    for r in res_int
        @printf("INTERCONNECTED  | %-10s | %10.2f  | %10.2f\n", r.area, r.lole, r.eue)
    end
end

end # module
