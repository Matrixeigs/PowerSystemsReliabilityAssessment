function [TestSystem, load_bus_indices] = dispaload(TestSystem)
%% ========================================================================
%  CONFIGURE DISPATCHABLE LOADS
%  ========================================================================
%  Purpose:
%  Modifies the MATPOWER case struct to model loads as "Virtual Generators".
%  
%  Why?
%  Standard OPF minimizes generation cost while meeting fixed load.
%  In reliability studies, we want to minimize *Load Shedding* when generation
%  is insufficient.
%
%  Methodology:
%  1. Convert each Load P_d into a Generator with P_g = -P_d.
%  2. Set P_min = -P_d (Max Load) and P_max = 0 (No Load).
%  3. Assign a high cost to these virtual generators.
%     - If the system is adequate, P_g stays at P_min (Full Load).
%     - If inadequate, P_g increases towards 0 (Load Shedding).
%     - The solver minimizes Cost = C * P_g. Since P_g is negative, 
%       this maximizes the absolute value of supplied load.
%
%  Inputs:
%  - TestSystem: Standard MATPOWER case struct
%
%  Outputs:
%  - TestSystem: Modified case with virtual generators
%  - load_bus_indices: Indices of buses where loads were located
%  ========================================================================

    %% 1. Identify Load Buses
    % Find buses with non-zero active power demand (Column 3)
    load_bus_indices = find(TestSystem.bus(:, 3) ~= 0);
    num_load_buses = length(load_bus_indices);
    
    % Store total original load for reference
    TestSystem.Pload = sum(TestSystem.bus(:, 3));

    %% 2. Configure Generator Costs
    % Real Generators: Set cost to 0 (or low) so they are used first.
    % Format: [model, startup, shutdown, n, c1, c0, ...]
    % Model 2 = Polynomial cost.
    
    num_real_gens = size(TestSystem.gen, 1);
    
    % Reset costs for real generators (Cost = 0)
    % This prioritizes using all available generation before shedding load.
    TestSystem.gencost = repmat([2 0 0 3 0 0 0], num_real_gens, 1);
    
    % Virtual Generators (Loads): Set cost to 1 (High)
    % Cost function: C(P) = 1 * P. 
    % Since P is negative (supply), minimizing 1*P tries to make P as negative as possible.
    % i.e., Maximize Supply.
    virtual_gen_cost = repmat([2 0 0 3 0 1 0], num_load_buses, 1);
    
    TestSystem.gencost = vertcat(TestSystem.gencost, virtual_gen_cost);

    %% 3. Create Virtual Generators
    % Initialize matrix for new generators
    % Size: [NumLoads x NumGenColumns]
    virtual_gens = zeros(num_load_buses, size(TestSystem.gen, 2));
    
    % Extract Load Data
    bus_ids = TestSystem.bus(load_bus_indices, 1);
    P_load  = TestSystem.bus(load_bus_indices, 3);
    Q_load  = TestSystem.bus(load_bus_indices, 4);
    
    % Map to Generator Columns (Standard MATPOWER format):
    % 1:BusID, 2:Pg, 3:Qg, 4:Qmax, 5:Qmin, 6:Vg, 7:mBase, 8:Status, 9:Pmax, 10:Pmin
    
    virtual_gens(:, 1)  = bus_ids;
    virtual_gens(:, 2)  = -P_load;          % Pg (Negative of Load)
    virtual_gens(:, 3)  = -Q_load;          % Qg
    virtual_gens(:, 4)  = 0;                % Qmax
    virtual_gens(:, 5)  = -Q_load;          % Qmin
    virtual_gens(:, 6)  = 1.0;              % Voltage Setpoint (dummy)
    virtual_gens(:, 7)  = TestSystem.baseMVA; % Base MVA
    virtual_gens(:, 8)  = 1;                % Status (1=On)
    virtual_gens(:, 9)  = 0;                % Pmax (0 = No Load Served / Full Shed)
    virtual_gens(:, 10) = -P_load;          % Pmin (Negative = Full Load Served)

    % Append to Generator Matrix
    TestSystem.gen = vertcat(TestSystem.gen, virtual_gens);

    %% 4. Clear Original Loads
    % Remove fixed loads from bus matrix to avoid double counting
    TestSystem.bus(load_bus_indices, 3:4) = 0;

end
