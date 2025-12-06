function [dns, nodal_dns] = mc_simulation(component_states, TestSystem, mpopt, numGenerators, numLines)

%% ========================================================================
%  MONTE CARLO SIMULATION: STATE EVALUATION
%  ========================================================================
%  Purpose:
%  Evaluates a specific system state (sampled in the Monte Carlo process)
%  using DC Optimal Power Flow (DC-OPF) to determine if load shedding
%  is required.
%
%  Inputs:
%  - component_states: Binary vector [1 x (Ng+Nl)] where 1=Failure, 0=Success
%  - TestSystem:       MATPOWER case structure
%  - mpopt:            MATPOWER options structure
%  - numGenerators:    Number of generators
%  - numLines:         Number of transmission lines
%
%  Outputs:
%  - dns:       Total Demand Not Supplied (MW) for this specific state
%  - nodal_dns: Vector [1 x NumBuses] of Demand Not Supplied per bus
%  ========================================================================

    % Disable warnings for singular matrices which often occur in failed states
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');

    %% 1. Map Sampled States to System Model
    % Extract generator and branch statuses from the input vector
    % Note: Input '1' means Failure, '0' means Normal
    
    % Generators: Column 8 in MATPOWER 'gen' matrix is status (1=ON, 0=OFF)
    gen_failures = component_states(:, 1:numGenerators);
    TestSystem.gen(1:numGenerators, 8) = ~gen_failures'; % If failure(1), status becomes OFF(0)

    % Branches: Column 11 in MATPOWER 'branch' matrix is status (1=ON, 0=OFF)
    branch_failures = component_states(:, numGenerators+1 : numGenerators+numLines);
    TestSystem.branch(1:numLines, 11) = ~branch_failures'; % If failure(1), status becomes OFF(0)

    %% 2. Run Optimal Power Flow (OPF)
    % Solve the optimization problem to minimize load shedding
    Result = runopf(TestSystem, mpopt);

    %% 3. Calculate Load Shedding (Demand Not Supplied)
    % Logic Explanation:
    % In nsqMain.m, loads are modeled as virtual generators with negative output (P_g < 0).
    % The cost function is defined as Cost = 1 * P_g.
    % - If load is fully served: P_g = -P_load  => Cost = -P_load
    % - If load is fully shed:   P_g = 0        => Cost = 0
    % 
    % Therefore: Result.f (Total Cost) = -P_served
    % DNS (Load Shed) = Total Load - P_served
    %                 = Total Load + Result.f
    
    dns = Result.f + TestSystem.load;

    % Filter out numerical noise (floating point errors)
    if dns < 0.1
        dns = 0;
    end
    
    %% 4. Calculate Nodal Curtailment (Weak Point Analysis Support)
    num_buses = size(TestSystem.bus, 1);
    nodal_dns = zeros(1, num_buses);
    
    if dns > 0
        % Virtual generators representing loads are appended after real generators
        total_gens_in_model = size(TestSystem.gen, 1);
        
        % Indices of virtual generators (loads)
        virt_gen_indices = (numGenerators + 1) : total_gens_in_model;
        
        if ~isempty(virt_gen_indices)
            % Extract data for virtual generators from OPF result
            % Col 1: Bus ID, Col 2: Pg (Output), Col 10: Pmin (Negative Load)
            virt_gens = Result.gen(virt_gen_indices, :);
            
            bus_ids = virt_gens(:, 1);
            P_gen   = virt_gens(:, 2);  % Actual generation (negative if supplying load)
            P_min   = virt_gens(:, 10); % Max capacity (negative of load demand)
            
            % Shed Amount = (Desired Supply) - (Actual Supply)
            % Desired Supply = -P_min (positive value)
            % Actual Supply  = -P_gen (positive value)
            % Shed = (-P_min) - (-P_gen) = P_gen - P_min
            
            shed_values = P_gen - P_min;
            
            % Map to nodal vector
            for k = 1:length(bus_ids)
                if shed_values(k) > 1e-3 % Threshold for noise
                    % Assuming bus IDs map directly to indices (1..24)
                    bus_idx = bus_ids(k);
                    if bus_idx <= num_buses
                        nodal_dns(bus_idx) = shed_values(k);
                    end
                end
            end
        end
    end
    
end


