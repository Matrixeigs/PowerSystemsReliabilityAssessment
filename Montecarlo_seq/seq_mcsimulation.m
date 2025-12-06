function [curtailment_mw, nodal_curtailment] = seq_mcsimulation(component_status, load_scale_factor, TestSystem, mpopt, load_bus_indices, numGenerators, numLines)
%% ========================================================================
%  HOURLY STATE EVALUATION (SEQUENTIAL)
%  ========================================================================
%  Purpose:
%  Evaluates the system for a specific hour in the chronological sequence.
%  1. Scales load based on the hourly factor.
%  2. Applies component outages based on the sampled state.
%  3. Runs DC-OPF to determine load curtailment.
%
%  Inputs:
%  - component_status:  Binary vector [Ng+Nl x 1] (1=Down, 0=Up)
%  - load_scale_factor: Scalar (0.0 to 1.0) representing % of Peak Load
%  - TestSystem:        MATPOWER case struct (with dispatchable loads)
%  - mpopt:             MATPOWER options
%  - load_bus_indices:  Indices of buses with loads (for virtual gen mapping)
%
%  Output:
%  - curtailment_mw:    Total Load Shedding (MW) for this hour
%  - nodal_curtailment: Vector [1 x NumBuses] of MW shed per bus
%  ========================================================================

    % Disable warnings for singular matrices which often occur in failed states
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');

    %% 1. Apply Chronological Load Scaling
    % The TestSystem passed here already has loads converted to virtual generators
    % (done in dispaload.m). These virtual generators are at indices > numGenerators.
    
    % Calculate indices for virtual generators representing loads
    virt_gen_start = numGenerators + 1;
    virt_gen_end   = numGenerators + length(load_bus_indices);
    
    % Scale the parameters of these virtual generators:
    % Col 2 (Pg), Col 3 (Qg), Col 5 (Qmin), Col 10 (Pmin/Pmax)
    % Note: Pmin is negative load. Scaling it scales the demand.
    TestSystem.gen(virt_gen_start:virt_gen_end, [2 3 5 10]) = ...
        TestSystem.gen(virt_gen_start:virt_gen_end, [2 3 5 10]) * load_scale_factor;
        
    % Update total system load tracker (for reference/calculation)
    TestSystem.Pload = TestSystem.Pload * load_scale_factor;

    %% 2. Apply Component Outages
    % Map binary status (1=Down) to MATPOWER status (0=Off)
    
    % Generators (Column 8)
    gen_failures = component_status(1:numGenerators);
    TestSystem.gen(1:numGenerators, 8) = ~gen_failures; 
    
    % Branches (Column 11)
    branch_failures = component_status(numGenerators+1 : numGenerators+numLines);
    TestSystem.branch(1:numLines, 11) = ~branch_failures;

    %% 3. Run DC-OPF
    % Solve optimization to minimize cost (which equals minimizing load shedding)
    Result = runopf(TestSystem, mpopt);
    
    %% 4. Calculate Load Curtailment
    % Objective Function (f) = Cost of Generation
    % Real Gen Cost = 0. Virtual Gen Cost = Penalty * MW_Supplied.
    % Since Virtual Gen MW is negative (load), MATPOWER minimizes:
    % Cost = Penalty * (Negative_Load_Supplied)
    %
    % DNS = Result.f + Total_Load
    
    curtailment_mw = Result.f + TestSystem.Pload;
    
    % Filter numerical noise
    if curtailment_mw < 0.1
        curtailment_mw = 0;
    end

    %% 5. Calculate Nodal Curtailment
    num_buses = size(TestSystem.bus, 1);
    nodal_curtailment = zeros(1, num_buses);
    
    if curtailment_mw > 0
        % Extract data for virtual generators (loads) from OPF result
        virt_gens = Result.gen(virt_gen_start:virt_gen_end, :);
        
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
                bus_idx = bus_ids(k);
                if bus_idx <= num_buses
                    nodal_curtailment(bus_idx) = shed_values(k);
                end
            end
        end
    end

end