% filepath: /Users/matrixeigs/Codes/PowerSystemsReliabilityAssessment/Montecarlo_seq/seq_mcsampling.m
function state_duration_matrix = seq_mcsampling(reliability_data, numGenerators, numLines, num_years, hours_per_year)
%% ========================================================================
%  SEQUENTIAL STATE SAMPLING (CHRONOLOGICAL)
%  ========================================================================
%  Purpose:
%  Generates a chronological sequence of component states (Up/Down) for
%  the specified duration.
%
%  Methodology:
%  Uses the "Next Event" method. For each component:
%  1. Start in Normal State (Up).
%  2. Sample Time-To-Failure (TTF) from Exp(lambda).
%  3. Sample Time-To-Repair (TTR) from Exp(mu).
%  4. Construct the time line by alternating these durations.
%
%  Inputs:
%  - reliability_data: Matrix [Ng+Nl x 2] containing [MTTF, MTTR]
%  - num_years:        Number of years to sample
%  - hours_per_year:   Hours in a simulation year (e.g., 8736 or 8760)
%
%  Outputs:
%  - state_duration_matrix: Sparse Matrix [(Ng+Nl) x TotalHours]
%                           1 = Component Down, 0 = Component Up
%  ========================================================================

    total_hours = num_years * hours_per_year;
    num_components = numGenerators + numLines;
    
    % Initialize sparse matrix to store failures (1s)
    % Using sparse is efficient because failures are rare events
    state_duration_matrix = sparse(num_components, total_hours);

    % Loop through each component to generate its history
    for i = 1 : num_components
        
        mttf = reliability_data(i, 1); % Mean Time To Failure
        mttr = reliability_data(i, 2); % Mean Time To Repair
        
        current_time = 0;
        is_up_state = true; % Assume components start in UP state
        
        % Generate timeline until we exceed total simulation hours
        while current_time < total_hours
            
            % Sample duration from Exponential Distribution
            % t = -Mean * ln(U), where U ~ Uniform(0,1)
            rand_val = rand(1);
            
            if is_up_state
                % Component is UP, sample Time To Failure
                duration = -mttf * log(rand_val);
                duration_int = round(duration); % Round to nearest hour
                
                % Advance time (no change to matrix, as 0 is default)
                current_time = current_time + duration_int;
            else
                % Component is DOWN, sample Time To Repair
                duration = -mttr * log(rand_val);
                duration_int = ceil(duration); % Ensure at least 1 hour repair
                
                % Mark these hours as Failure (1) in the matrix
                start_idx = round(current_time) + 1;
                end_idx = min(start_idx + duration_int - 1, total_hours);
                
                if start_idx <= total_hours
                    state_duration_matrix(i, start_idx:end_idx) = 1;
                end
                
                current_time = current_time + duration_int;
            end
            
            % Toggle state for next iteration
            is_up_state = ~is_up_state;
        end
    end
end












