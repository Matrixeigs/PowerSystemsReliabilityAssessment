%% ========================================================================
%  SEQUENTIAL MONTE CARLO SIMULATION (HL2)
%  ========================================================================
%  Educational Implementation for PhD Students
%
%  PURPOSE:
%  Performs reliability assessment using Sequential Monte Carlo Simulation (SMCS).
%  Unlike Non-Sequential MCS (State Sampling), SMCS simulates the system
%  chronologically (hour by hour) over many years.
%
%  KEY CONCEPTS:
%  1. Chronological Load Model: Load varies hourly based on IEEE RTS profile.
%  2. State Duration Sampling: Component states (Up/Down) are sampled based
%     on probability distributions of Time-To-Failure (TTF) and Time-To-Repair (TTR).
%  3. Frequency & Duration Indices: SMCS allows calculating frequency (occ/yr)
%     and duration (hours/disturbance) indices, which State Sampling cannot easily do.
%  ========================================================================

clc; close all; clear;
tic;

% Disable warnings globally for cleaner output
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

%% ========================================================================
%  SECTION 1: INITIALIZATION
%  ========================================================================
fprintf('Initializing Sequential Simulation...\n');

% Load System Data
TestSystem = loadcase('case24_ieee_rts');
numBuses = size(TestSystem.bus, 1);
numGenerators = size(TestSystem.gen, 1);
numLines = size(TestSystem.branch, 1);

% Simulation Parameters
HOURS_PER_YEAR = 8736;      % IEEE RTS uses 52 weeks * 7 days * 24 hours
MAX_SIM_YEARS = 4000;       % Maximum simulation horizon
COV_THRESHOLD = 0.05;       % Convergence criterion (Coefficient of Variation)
CURTAIL_THRESHOLD = 0.01;   % MW threshold to record a curtailment event

% Initialize Storage for Indices (Per Year)
results_year.plc = zeros(MAX_SIM_YEARS, 1);  % Probability of Load Curtailment
results_year.nlc = zeros(MAX_SIM_YEARS, 1);  % Number of Load Curtailments (Freq)
results_year.dlc = zeros(MAX_SIM_YEARS, 1);  % Duration of Load Curtailments
results_year.dns = zeros(MAX_SIM_YEARS, 1);  % Demand Not Supplied (MW avg)
results_year.ens = zeros(MAX_SIM_YEARS, 1);  % Energy Not Supplied (MWh)
results_year.bpii = zeros(MAX_SIM_YEARS, 1); % Bulk Power Interruption Index

% Initialize Storage for Cumulative Indices (Running Average)
results_cum.eens = zeros(MAX_SIM_YEARS, 1);
results_cum.cov = zeros(MAX_SIM_YEARS, 1);

% Initialize Accumulators for Post-Processing
accum_nodal_eens = zeros(1, numBuses); % Accumulate MWh shed per bus
accum_comp_fail_during_loss = zeros(numGenerators + numLines, 1); % For weak point detection
total_loss_hours = 0; % Total hours where load shedding occurred

%% ========================================================================
%  SECTION 2: CHRONOLOGICAL MODELING
%  ========================================================================
fprintf('Generating Chronological Models...\n');

% 1. Generate Hourly Load Curve (8736 hours)
% Returns: Active Load, Reactive Load, and Scaling Factors per hour
[busPd, busQd, load_scale_factors] = anloducurve(HOURS_PER_YEAR); 

% 2. Prepare Component Reliability Data (MTTF, MTTR)
% Used to generate chronological Up/Down sequences
reliability_data = seqmeantime();

% 3. Configure Dispatchable Loads for OPF
% Modifies the MATPOWER case to treat loads as negative generators with high cost
[TestSystem, load_bus_indices] = dispaload(TestSystem);

% 4. OPF Options
mpopt = mpoption('PF_DC',1, 'VERBOSE',0, 'OUT_ALL',0, 'OPF_ALG_DC',200, 'OPF_FLOW_LIM',1);

%% ========================================================================
%  SECTION 3: SEQUENTIAL MONTE CARLO LOOP
%  ========================================================================
fprintf('Starting Simulation Loop...\n');

for iYear = 1 : MAX_SIM_YEARS
    
    %% --- Step 3.1: Sequential State Sampling ---
    % Generate a time series (matrix) of component states for the entire year.
    % Matrix size: (Ng+Nl) x HOURS_PER_YEAR
    % 1 = Component Down, 0 = Component Up
    system_state_chronological = seq_mcsampling(reliability_data, numGenerators, numLines, 1, HOURS_PER_YEAR);

    %% --- Step 3.2: Identify Contingency Hours ---
    % Optimization: We only need to run OPF if there is a component failure.
    % (Assuming N-0 system is secure).
    % Find columns (hours) where at least one component is down (sum > 0).
    contingency_hours = find(sum(system_state_chronological, 1) > 0);
    
    % Extract states for these specific hours
    contingency_states = system_state_chronological(:, contingency_hours);
    
    %% --- Step 3.3: Chronological Evaluation (OPF) ---
    num_contingencies = length(contingency_hours);
    hourly_curtailment = zeros(1, num_contingencies);
    hourly_nodal_curtailment = zeros(num_contingencies, numBuses);
    
    if num_contingencies > 0
        fprintf('  Year %d: Found %d contingency hours. Starting parallel evaluation...\n', iYear, num_contingencies);
    end
    
    % Parallel evaluation of contingency hours
    parfor j = 1 : num_contingencies
        hour_idx = contingency_hours(j);
        current_load_scale = load_scale_factors(hour_idx);
        current_state = contingency_states(:, j);
        
        try
            % Run OPF for this specific hour's load level and component status
            [curtailment, nodal_curt] = seq_mcsimulation(current_state, current_load_scale, TestSystem, mpopt, load_bus_indices, numGenerators, numLines);
            hourly_curtailment(j) = curtailment;
            hourly_nodal_curtailment(j, :) = nodal_curt;
        catch
            % Fallback for solver failures
            hourly_curtailment(j) = 0; 
            hourly_nodal_curtailment(j, :) = zeros(1, numBuses);
        end
        
        % Progress monitoring: Print every ~20% of progress
        % Note: In parfor, execution order is not guaranteed, so indices may appear out of order.
        if mod(j, max(1, round(num_contingencies/5))) == 0
             fprintf('    -> [Year %d] Worker processed contingency index %d / %d\n', iYear, j, num_contingencies);
        end
    end
    
    % Create a full-year curtailment vector (mostly zeros)
    year_curtailment_profile = zeros(1, HOURS_PER_YEAR);
    year_curtailment_profile(contingency_hours) = hourly_curtailment;
    
    % Binary flag profile (1 if curtailment > 0)
    year_curtailment_flag = double(year_curtailment_profile > CURTAIL_THRESHOLD);

    %% --- Step 3.4: Update Accumulators for Post-Processing ---
    % Identify indices within the contingency set that resulted in loss
    loss_indices_in_contingency = find(hourly_curtailment > CURTAIL_THRESHOLD);
    
    if ~isempty(loss_indices_in_contingency)
        % 1. Accumulate Nodal EENS (MWh)
        % Sum nodal curtailment for all loss hours in this year
        accum_nodal_eens = accum_nodal_eens + sum(hourly_nodal_curtailment(loss_indices_in_contingency, :), 1);
        
        % 2. Accumulate Component Failures during Loss (Weak Point Detection)
        % Get states for hours where loss occurred
        states_during_loss = contingency_states(:, loss_indices_in_contingency);
        % Sum failures per component
        accum_comp_fail_during_loss = accum_comp_fail_during_loss + sum(states_during_loss, 2);
        
        total_loss_hours = total_loss_hours + length(loss_indices_in_contingency);
    end

    %% --- Step 3.5: Calculate Annual Indices ---
    % PLC: Probability (Fraction of time load is curtailed)
    results_year.plc(iYear) = sum(year_curtailment_flag) / HOURS_PER_YEAR;
    
    % NLC: Number of Load Curtailments (Frequency)
    % Uses a helper to count distinct events in the time series
    results_year.nlc(iYear) = calnlc(year_curtailment_flag);
    
    % DLC: Duration of Load Curtailment (Total hours)
    results_year.dlc(iYear) = sum(year_curtailment_flag);
    
    % EENS: Energy Not Supplied (MWh)
    % Sum of MW curtailed in each hour * 1 hour
    results_year.ens(iYear) = sum(year_curtailment_profile);
    
    % EDNS: Average Demand Not Supplied (MW)
    results_year.dns(iYear) = results_year.ens(iYear) / HOURS_PER_YEAR;

    %% --- Step 3.6: Update Expected Indices & Check Convergence ---
    % Calculate running mean of EENS
    results_cum.eens(iYear) = mean(results_year.ens(1:iYear));
    
    % Calculate Coefficient of Variation (Standard Error / Mean)
    if iYear > 1
        std_dev = std(results_year.ens(1:iYear));
        results_cum.cov(iYear) = std_dev / (results_cum.eens(iYear) * sqrt(iYear));
        
        % Display progress
        if mod(iYear, 10) == 0
            fprintf('Year %4d | EENS: %.4f MWh/yr | CoV: %.4f\n', ...
                iYear, results_cum.eens(iYear), results_cum.cov(iYear));
        end
        
        % Check Convergence
        if results_cum.cov(iYear) < COV_THRESHOLD && results_cum.cov(iYear) > 0
            fprintf('Convergence Reached at Year %d!\n', iYear);
            break;
        end
    end
end

% Trim unused pre-allocated space
final_year = iYear;
results_cum.eens = results_cum.eens(1:final_year);
results_cum.cov = results_cum.cov(1:final_year);
% ... (trim others similarly if needed)

%% ========================================================================
%  SECTION 4: FINAL RESULTS & POST-PROCESSING
%  ========================================================================
fprintf('\n--- SIMULATION COMPLETE ---\n');
fprintf('EENS (Expected Energy Not Supplied): %.4f MWh/yr\n', results_cum.eens(end));
fprintf('LOLE (Loss of Load Expectation):     %.4f hr/yr\n', mean(results_year.dlc(1:final_year)));
fprintf('LOLF (Loss of Load Frequency):       %.4f occ/yr\n', mean(results_year.nlc(1:final_year)));

%% --- 4.1 Nodal Analysis ---
fprintf('\n--- NODAL RELIABILITY INDICES ---\n');
% Calculate Average Annual Nodal EENS
nodal_eens_avg = accum_nodal_eens / final_year; % MWh/yr

% Display top 5 worst buses
[sorted_eens, sorted_idx] = sort(nodal_eens_avg, 'descend');
fprintf('Top 5 Buses by EENS (MWh/yr):\n');
for k = 1:min(5, length(sorted_eens))
    if sorted_eens(k) > 0
        fprintf('  Bus %2d: %.4f MWh/yr\n', sorted_idx(k), sorted_eens(k));
    end
end

%% --- 4.2 Weak Point Detection ---
fprintf('\n--- WEAK POINT DETECTION ---\n');
% Metric: Probability(Component Down | System Failure)
if total_loss_hours > 0
    comp_importance = accum_comp_fail_during_loss / total_loss_hours;
    
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

%% --- 4.3 Export Results ---
fprintf('\n--- EXPORTING DATA ---\n');
% Export Nodal Results to CSV
bus_ids = (1:numBuses)';
nodal_table = table(bus_ids, nodal_eens_avg', 'VariableNames', {'BusID', 'EENS_MWh_yr'});
writetable(nodal_table, 'seq_nodal_results.csv');
fprintf('Saved nodal results to "seq_nodal_results.csv"\n');

% Save full workspace
save('seq_reliability_results.mat', 'results_year', 'results_cum', ...
     'nodal_eens_avg', 'comp_importance');
fprintf('Saved full results to "seq_reliability_results.mat"\n');

fprintf('========================================\n\n');

%% ========================================================================
%  SECTION 5: VISUALIZATION
%  ========================================================================
fprintf('Generating plots...\n');

figure('Name', 'Sequential MCS Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

% Subplot 1: Convergence of EENS
subplot(2, 3, 1);
plot(1:final_year, results_cum.eens, 'b-', 'LineWidth', 1.5);
title('Convergence of EENS'); xlabel('Simulation Year'); ylabel('MWh/yr'); grid on;

% Subplot 2: Convergence of CoV
subplot(2, 3, 2);
plot(1:final_year, results_cum.cov, 'r-', 'LineWidth', 1.5);
yline(COV_THRESHOLD, 'k--', 'Target');
title('Convergence of CoV'); xlabel('Simulation Year'); grid on;

% Subplot 3: Annual Variability (Histogram of EENS)
subplot(2, 3, 3);
histogram(results_year.ens(1:final_year), 20, 'FaceColor', [0.4 0.4 0.4]);
title('Distribution of Annual EENS'); xlabel('MWh/yr'); ylabel('Frequency'); grid on;

% Subplot 4: Nodal EENS
subplot(2, 3, 4);
bar(nodal_eens_avg, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('Bus ID'); ylabel('EENS (MWh/yr)');
title('Nodal Reliability'); grid on; xlim([0 numBuses+1]);

% Subplot 5: Critical Components
subplot(2, 3, [5 6]);
num_top = min(15, length(comp_importance));
[top_imp, top_idx] = maxk(comp_importance, num_top);
bar(top_imp * 100, 'FaceColor', [0.8, 0.3, 0.3]);
xticks(1:num_top);
labels = cell(num_top, 1);
for k = 1:num_top
    idx = top_idx(k);
    if idx <= numGenerators
        labels{k} = sprintf('G%d', idx);
    else
        labels{k} = sprintf('L%d', idx - numGenerators);
    end
end
xticklabels(labels); xtickangle(45);
ylabel('Probability (%)'); title(['Top ', num2str(num_top), ' Critical Components']); grid on;

fprintf('Visualization complete.\n');

toc;
