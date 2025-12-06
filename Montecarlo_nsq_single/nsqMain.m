% filepath: /Users/matrixeigs/Codes/PowerSystemsReliabilityAssessment/Montecarlo_nsq_single/nsqMain.m
%% ========================================================================
%  MONTE CARLO NON-SEQUENTIAL SIMULATION FOR POWER SYSTEM RELIABILITY
%  ========================================================================
%  Educational Implementation for PhD Students
%  
%  PURPOSE:
%  This code performs Non-Sequential Monte Carlo (NSQ) simulation to assess
%  bulk power system reliability using IEEE RTS-24 bus test system.
%
%  RELIABILITY INDICES COMPUTED:
%  - ENLC:  Expected Number of Load Curtailments (occ./yr)
%  - EDLC:  Expected Duration of Load Curtailments (hr/yr)
%  - PLC:   Probability of Load Curtailments
%  - EDNS:  Expected Demand Not Supplied (MW)
%  - EENS:  Expected Energy Not Supplied (MWh/yr)
%  - BPII:  Bulk Power Interruption Index (MW/MW-yr)
%  - BPECI: Bulk Power Energy Curtailment Index (MWh/MW-yr)
%  - SI:    Severity Index (system minutes/yr)
%
%  METHODOLOGY:
%  1. Monte Carlo sampling of component states (generators & lines)
%  2. DC Optimal Power Flow for each sampled state
%  3. Statistical convergence using coefficient of variation (beta)
%  ========================================================================

clear; close all; clc;
tic; % Start total execution timer

% Disable warnings globally for cleaner output
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

%% ========================================================================
%  SECTION 1: SYSTEM INITIALIZATION
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 1: Loading Test System\n');
fprintf('========================================\n');

% Load IEEE RTS-24 bus test system
TestSystem = loadcase('case24_ieee_rts');
numBuses = size(TestSystem.bus, 1);
numGenerators = size(TestSystem.gen, 1);
numLines = size(TestSystem.branch, 1);

fprintf('System loaded: %d buses, %d generators, %d lines\n', ...
        numBuses, numGenerators, numLines);

%% ========================================================================
%  SECTION 2: SIMULATION PARAMETERS
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 2: Setting Simulation Parameters\n');
fprintf('========================================\n');

% Convergence criterion: Coefficient of Variation (beta)
% beta = std(EDNS) / mean(EDNS)
% Simulation stops when beta < beta_limit
beta_limit = 0.0017;          % Target accuracy (0.17%)
max_iterations = 100000;      % Maximum number of samples
samples_per_batch = 100;      % Samples processed in each iteration

fprintf('Beta convergence limit: %.4f\n', beta_limit);
fprintf('Maximum iterations: %d\n', max_iterations);
fprintf('Samples per batch: %d\n', samples_per_batch);

% Initialize iteration counters and indices
current_iteration = 0;
accumulated_edns = 0;         % Expected Demand Not Supplied
accumulated_lole = 0;         % Loss of Load Expectation
current_beta = inf;           % Coefficient of variation
database_row_count = 0;       % Tracks rows in state database

%% ========================================================================
%  SECTION 3: PRE-ALLOCATE DATA STRUCTURES
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 3: Allocating Memory\n');
fprintf('========================================\n');

% Educational Note: Pre-allocation improves performance in MATLAB
% Each row represents a unique system state with:
% Columns 1 to numGenerators: Generator status (1=up, 0=down)
% Columns numGenerators+1 to numGenerators+numLines: Line status
% Column numGenerators+numLines+1: Occurrence count
% Column numGenerators+numLines+2: Load curtailment (MW)
% Column numGenerators+numLines+3: Failure indicator (1=failure, 0=success)
% Columns ...+4 to ...+3+numBuses: Nodal Load Curtailment (MW) [NEW]

col_count = numGenerators + numLines + 1;
col_dns   = numGenerators + numLines + 2;
col_flag  = numGenerators + numLines + 3;
col_nodal_start = numGenerators + numLines + 4;
col_nodal_end   = numGenerators + numLines + 3 + numBuses;

total_cols = col_nodal_end;

state_database = sparse(max_iterations, total_cols);

% Convergence tracking arrays
num_checkpoints = max_iterations / samples_per_batch;
beta_history = zeros(1, num_checkpoints);
edns_history = zeros(1, num_checkpoints);
lole_history = zeros(1, num_checkpoints);
plc_history = zeros(1, num_checkpoints);

fprintf('Memory allocated for %d states\n', max_iterations);

%% ========================================================================
%  SECTION 4: LOAD MODELING (DISPATCHABLE LOAD METHOD)
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 4: Modeling Loads as Negative Generators\n');
fprintf('========================================\n');

% Educational Note: In DC-OPF, we model load curtailment as "negative generation"
% This allows the optimizer to shed load when generation is insufficient

% Identify buses with loads
load_buses = find(TestSystem.bus(:, 3) ~= 0);
num_load_buses = size(load_buses, 1);

% Store total system load
TestSystem.load = sum(TestSystem.bus(:, 3));

% Set all real generator costs to zero (we only penalize load shedding)
TestSystem.gencost = repmat([2 0 0 3 0 0 0], numGenerators, 1);

% Create virtual generators at load buses with high cost
% Cost = 1 ensures load shedding is minimized in OPF
TestSystem.gencost = vertcat(TestSystem.gencost, ...
                             repmat([2 0 0 3 0 1 0], num_load_buses, 1));

% Configure virtual generator parameters
virtual_gen_template = TestSystem.gen(1:num_load_buses, :);
virtual_gen_template(:, 1:10) = [ ...
    TestSystem.bus(load_buses, 1), ...                    % Bus number
    -TestSystem.bus(load_buses, 3), ...                   % Pg (negative = supply)
    -TestSystem.bus(load_buses, 4), ...                   % Qg
    zeros(num_load_buses, 1), ...                         % Qmax
    -TestSystem.bus(load_buses, 4), ...                   % Qmin
    zeros(num_load_buses, 1), ...                         % Vg
    TestSystem.baseMVA * ones(num_load_buses, 1), ...     % mBase
    ones(num_load_buses, 1), ...                          % status
    zeros(num_load_buses, 1), ...                         % Pmax
    -TestSystem.bus(load_buses, 3) ...                    % Pmin (negative)
];

TestSystem.gen = [TestSystem.gen; virtual_gen_template];

% Remove original loads from buses (now represented as generators)
TestSystem.bus(load_buses, 3:4) = 0;

fprintf('Total system load: %.2f MW\n', TestSystem.load);
fprintf('Number of load points: %d\n', num_load_buses);

%% ========================================================================
%  SECTION 5: COMPONENT FAILURE PROBABILITIES
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 5: Loading Failure Probabilities\n');
fprintf('========================================\n');

% Load component failure probabilities from external function
% Returns vector: [generator probabilities; line probabilities]
component_failure_probs = failprob();

fprintf('Failure probabilities loaded for %d components\n', ...
        length(component_failure_probs));

%% ========================================================================
%  SECTION 6: OPTIMAL POWER FLOW SETTINGS
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 6: Configuring OPF Solver\n');
fprintf('========================================\n');

% MATPOWER options for DC-OPF
% PF_DC=1: Use DC power flow
% VERBOSE=0: Suppress solver output
% OUT_ALL=0: Suppress result output
% OPF_ALG_DC=200: Use MIPS solver for DC-OPF
% OPF_FLOW_LIM=1: Enforce branch flow limits
mpopt = mpoption('PF_DC', 1, 'VERBOSE', 0, 'OUT_ALL', 0, ...
                 'OPF_ALG_DC', 200, 'OPF_FLOW_LIM', 1);

% Run baseline OPF to verify system setup
baseline_result = runopf(TestSystem, mpopt);

fprintf('OPF solver configured (DC approximation)\n');

if baseline_result.success
    status_str = 'SUCCESS';
else
    status_str = 'FAILED';
end
fprintf('Baseline system check: %s\n', status_str);

%% ========================================================================
%  SECTION 7: MAIN MONTE CARLO LOOP
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 7: Starting Monte Carlo Simulation\n');
fprintf('========================================\n');
fprintf('Progress will be displayed every %d samples\n\n', samples_per_batch);

while current_beta > beta_limit && current_iteration < max_iterations
    
    %% --- Step 7.1: Monte Carlo Sampling ---
    % Generate random component states based on failure probabilities
    sampled_states = mc_sampling(component_failure_probs, samples_per_batch, ...
                                 numGenerators, numLines);
    
    % Initialize columns for: count, load_curtailed, failure_flag
    % Note: mc_sampling returns [samples x (Ng+Nl)]
    % We need to append columns for Count, DNS, Flag, and Nodal Results
    
    % Remove duplicate states and accumulate their counts
    [sampled_states, ~, duplicate_indices] = unique(sampled_states, 'rows', 'stable');
    
    % Calculate counts
    counts = zeros(size(sampled_states, 1), 1);
    for i = 1:size(sampled_states, 1)
        counts(i) = sum(duplicate_indices == i);
    end
    
    % Append Count column
    sampled_states = [sampled_states, counts]; 
    
    %% --- Step 7.2: Check Against Existing Database ---
    if current_iteration > 0
        % Find states already in database (compare only component statuses)
        [~, existing_rows, new_rows] = intersect(...
            state_database(1:database_row_count, 1:numGenerators+numLines), ...
            sampled_states(:, 1:numGenerators+numLines), 'rows', 'stable');
        
        % Accumulate counts for existing states
        state_database(existing_rows, col_count) = ...
            state_database(existing_rows, col_count) + ...
            sampled_states(new_rows, end); % The last column is now Count
        
        % Remove already-processed states
        sampled_states(new_rows, :) = [];
    end
    
    %% --- Step 7.3: Simulate New States ---
    % Educational Note: This is the computationally intensive part
    % For each new state, we solve DC-OPF to find minimum load curtailment
    
    num_new_states = size(sampled_states, 1);
    
    if num_new_states > 0
        % Pre-allocate results: [TotalDNS, NodalDNS_1...N]
        sim_results_dns = zeros(num_new_states, 1 + numBuses);
        
        parfor i = 1:num_new_states
            % Returns load curtailed (MW) and Nodal breakdown
            [d, nd] = mc_simulation(...
                sampled_states(i, 1:numGenerators+numLines), ...
                TestSystem, mpopt, numGenerators, numLines);
            sim_results_dns(i, :) = [d, nd];
        end
        
        % Construct the full rows to append to database
        % sampled_states currently has [Components, Count]
        % We need [Components, Count, DNS, Flag, NodalDNS...]
        
        total_dns = sim_results_dns(:, 1);
        fail_flag = double(total_dns > 1e-4);
        nodal_res = sim_results_dns(:, 2:end);
        
        new_database_entries = [sampled_states, total_dns, fail_flag, nodal_res];
        
        % Add new states to database
        state_database(database_row_count+1:database_row_count+num_new_states, :) = ...
            new_database_entries;
        database_row_count = database_row_count + num_new_states;
    end
    
    %% --- Step 7.4: Compute Reliability Indices ---
    current_samples = current_iteration + samples_per_batch;
    active_data = state_database(1:database_row_count, :);
    
    % EDNS: Expected Demand Not Supplied (MW)
    accumulated_edns = full(sum(active_data(:, col_count) .* ...
                          active_data(:, col_dns)) / current_samples);
    
    % LOLE: Loss of Load Expectation (hours/year)
    accumulated_lole = full(sum(active_data(:, col_count) .* ...
                          active_data(:, col_flag)) / ...
                          current_samples * 8760);
    
    % PLC: Probability of Load Curtailment
    plc = full(sum(active_data(:, col_count) .* ...
             active_data(:, col_flag)) / current_samples);
    
    % Beta: Coefficient of Variation (convergence measure)
    current_beta = full(sqrt(sum(active_data(:, col_count) .* ...
                           (active_data(:, col_dns) - accumulated_edns).^2)) / ...
                           current_samples / accumulated_edns);
    
    %% --- Step 7.5: Record Progress ---
    checkpoint_index = current_samples / samples_per_batch;
    beta_history(checkpoint_index) = current_beta;
    edns_history(checkpoint_index) = accumulated_edns;
    lole_history(checkpoint_index) = accumulated_lole;
    plc_history(checkpoint_index) = plc;
    
    % Update iteration counter
    current_iteration = current_samples;
    
    % Display progress every 1000 samples
    if mod(current_iteration, 1000) == 0
        fprintf('Iteration %6d: Beta = %.6f, EDNS = %.4f MW, LOLE = %.4f hr/yr\n', ...
                current_iteration, current_beta, accumulated_edns, accumulated_lole);
    end
end

elapsed_time = toc;

%% ========================================================================
%  SECTION 8: FINAL RESULTS & POST-PROCESSING
%  ========================================================================
fprintf('\n========================================\n');
fprintf('SECTION 8: SIMULATION RESULTS\n');
fprintf('========================================\n');
fprintf('Total simulation time: %.2f seconds\n', elapsed_time);
fprintf('Total iterations: %d\n', current_iteration);
fprintf('Unique states evaluated: %d\n', database_row_count);

if current_beta <= beta_limit
    conv_str = 'YES';
else
    conv_str = 'NO';
end
fprintf('Convergence achieved: %s\n', conv_str);

fprintf('\n--- RELIABILITY INDICES ---\n');
fprintf('EDNS (Expected Demand Not Supplied): %.4f MW\n', accumulated_edns);
fprintf('LOLE (Loss of Load Expectation): %.4f hours/year\n', accumulated_lole);
fprintf('PLC (Probability of Load Curtailment): %.6f\n', plc);
fprintf('Beta (Coefficient of Variation): %.6f\n', current_beta);

%% --- 8.1 Nodal Analysis ---
fprintf('\n--- NODAL RELIABILITY INDICES ---\n');
% Calculate EENS per bus: Sum(Count * NodalDNS) / TotalSamples
nodal_eens = full(sum(active_data(:, col_count) .* ...
                 active_data(:, col_nodal_start:col_nodal_end), 1) / current_samples);

% Display top 5 worst buses
[sorted_eens, sorted_idx] = sort(nodal_eens, 'descend');
fprintf('Top 5 Buses by EENS (MWh/yr):\n');
for k = 1:min(5, length(sorted_eens))
    if sorted_eens(k) > 0
        fprintf('  Bus %2d: %.4f MWh/yr\n', sorted_idx(k), sorted_eens(k) * 8760);
    end
end

%% --- 8.2 Weak Point Detection ---
fprintf('\n--- WEAK POINT DETECTION ---\n');
% Identify components that contribute most to system failures
% Metric: Probability(Component Down | System Failure)

% Filter rows where system failure occurred
failure_indices = find(active_data(:, col_flag) == 1);
if ~isempty(failure_indices)
    failure_data = active_data(failure_indices, :);
    failure_weights = failure_data(:, col_count); % Number of occurrences
    total_failure_events = sum(failure_weights);
    
    % Component states during failure (1=Down)
    comp_states_in_failure = failure_data(:, 1:(numGenerators+numLines));
    
    % Weighted sum of component failures during system failure
    comp_importance = full((comp_states_in_failure' * failure_weights) / total_failure_events);
    
    [sorted_imp, sorted_comp_idx] = sort(comp_importance, 'descend');
    
    fprintf('Top 5 Critical Components (Prob. Down given System Failure):\n');
    for k = 1:min(5, length(sorted_imp))
        c_idx = sorted_comp_idx(k);
        if c_idx <= numGenerators
            type = 'Gen'; id = c_idx;
        else
            type = 'Line'; id = c_idx - numGenerators;
        end
        fprintf('  %s %2d: %.2f%%\n', type, id, sorted_imp(k) * 100);
    end
else
    fprintf('No failure events recorded to analyze weak points.\n');
    comp_importance = zeros(numGenerators+numLines, 1);
end

%% --- 8.3 Export Results ---
fprintf('\n--- EXPORTING DATA ---\n');
% Export Nodal Results to CSV
bus_ids = (1:numBuses)';
nodal_table = table(bus_ids, nodal_eens' * 8760, 'VariableNames', {'BusID', 'EENS_MWh_yr'});
writetable(nodal_table, 'nodal_results.csv');
fprintf('Saved nodal results to "nodal_results.csv"\n');

% Save full workspace
save('reliability_results.mat', 'accumulated_edns', 'accumulated_lole', ...
     'nodal_eens', 'comp_importance', 'beta_history', 'edns_history');
fprintf('Saved full results to "reliability_results.mat"\n');

fprintf('========================================\n\n');

%% ========================================================================
%  SECTION 9: VISUALIZATION
%  ========================================================================
fprintf('Generating plots...\n');

actual_checkpoints = current_iteration / samples_per_batch;

% Figure 1: Convergence (Existing)
figure('Name', 'Convergence Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 1);
plot(1:actual_checkpoints, beta_history(1:actual_checkpoints), 'b-', 'LineWidth', 1.5);
hold on; yline(beta_limit, 'r--', 'LineWidth', 1.5);
title('Convergence of \beta'); xlabel('Iter (x100)'); grid on;

subplot(2, 2, 2);
plot(1:actual_checkpoints, edns_history(1:actual_checkpoints), 'g-', 'LineWidth', 1.5);
title('EDNS (MW)'); xlabel('Iter (x100)'); grid on;

subplot(2, 2, 3);
plot(1:actual_checkpoints, lole_history(1:actual_checkpoints), 'm-', 'LineWidth', 1.5);
title('LOLE (hr/yr)'); xlabel('Iter (x100)'); grid on;

subplot(2, 2, 4);
plot(1:actual_checkpoints, plc_history(1:actual_checkpoints), 'c-', 'LineWidth', 1.5);
title('Prob. Load Curtailment'); xlabel('Iter (x100)'); grid on;

% Figure 2: Nodal & Weak Point Analysis (New)
figure('Name', 'System Weak Points', 'NumberTitle', 'off', 'Position', [150, 150, 1000, 500]);

% Subplot 1: Nodal EENS
subplot(1, 2, 1);
bar(nodal_eens * 8760, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('Bus ID'); ylabel('EENS (MWh/yr)');
title('Nodal Reliability (EENS)');
grid on; xlim([0 numBuses+1]);

% Subplot 2: Component Importance
subplot(1, 2, 2);
% Plot only top 15 critical components for clarity
num_top = min(15, length(comp_importance));
[top_imp, top_idx] = maxk(comp_importance, num_top);
bar(top_imp * 100, 'FaceColor', [0.8, 0.3, 0.3]);
xticks(1:num_top);
% Create labels
labels = cell(num_top, 1);
for k = 1:num_top
    idx = top_idx(k);
    if idx <= numGenerators
        labels{k} = sprintf('G%d', idx);
    else
        labels{k} = sprintf('L%d', idx - numGenerators);
    end
end
xticklabels(labels);
xtickangle(45);
ylabel('Probability (%)');
title(['Top ', num2str(num_top), ' Critical Components']);
grid on;

fprintf('Visualization complete.\n');

%% ========================================================================
%  END OF SIMULATION
%  ========================================================================