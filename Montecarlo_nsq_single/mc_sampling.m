% filepath: /Users/matrixeigs/Codes/PowerSystemsReliabilityAssessment/Montecarlo_nsq_single/mc_sampling.m
function eqstatus = mc_sampling(failure_probabilities, num_samples, numGenerators, numLines)
%% ========================================================================
%  MONTE CARLO SAMPLING ENGINE
%  ========================================================================
%  Purpose:
%  Generates random system states based on component failure probabilities.
%  Uses the Inverse Transform Method for binary states (Bernoulli trials).
%
%  Inputs:
%  - failure_probabilities: Vector of unavailability (U) for all components
%  - num_samples:           Number of Monte Carlo years/states to generate
%  - numGenerators:         Number of generators
%  - numLines:              Number of transmission lines
%
%  Outputs:
%  - eqstatus: Sparse matrix [num_samples x (Ng+Nl)]
%              1 = Component Failed
%              0 = Component Available
%  ========================================================================

    %% 1. Generate Random Numbers
    % Create a matrix of uniform random numbers U ~ [0, 1]
    random_numbers = rand(num_samples, numGenerators + numLines);
    
    %% 2. Determine Component Status
    % Compare random numbers against failure probabilities (Unavailability)
    % If rand < Probability_of_Failure, then state is Failure (1)
    % This creates a Bernoulli distribution for each component
    
    % Expand probability vector to match matrix dimensions
    prob_matrix = repmat(failure_probabilities, num_samples, 1);
    
    % Perform comparison
    eqstatus = random_numbers < prob_matrix;
    
    %% 3. Apply Special Constraints
    % The synchronous compensator at Bus 14 (Index 15) is assumed perfectly reliable
    % for this study, or its failure is neglected.
    SYNC_COMP_INDEX = 15;
    eqstatus(:, SYNC_COMP_INDEX) = 0; 

    %% 4. Convert to Sparse Matrix
    % Use sparse storage to save memory, as failures are rare events (mostly 0s)
    eqstatus = sparse(double(eqstatus));

end


