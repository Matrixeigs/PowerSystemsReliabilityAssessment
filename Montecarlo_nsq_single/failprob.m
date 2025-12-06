function total_unavailability = failprob()
%% ========================================================================
%  COMPONENT UNAVAILABILITY CALCULATION
%  ========================================================================
%  Purpose:
%  Calculates the probability that a component is in a failed state 
%  (Unavailability, U) based on its failure rate and repair time.
%
%  Formulas:
%  1. Unavailability (U) = MTTR / (MTTF + MTTR)
%     where:
%       MTTF = Mean Time To Failure = 1 / lambda
%       MTTR = Mean Time To Repair  = 1 / mu
%
%  2. Alternatively: U = lambda / (lambda + mu)
%  ========================================================================

    % Load raw data
    data = case24_failrate();

    %% 1. Generator Unavailability
    % Using Formula 1: U = MTTR / (MTTF + MTTR)
    prob_gen = data.genmttr ./ (data.genmttf + data.genmttr);

    % Note: Planned maintenance is currently commented out in the original logic.
    % If included, it would add a scheduled outage rate component.

    %% 2. Branch Unavailability
    % First, calculate Repair Rate (mu) from Duration (r)
    % mu = 8760 hours/year / Duration (hours)
    branch_mu = 8760 ./ data.brdur;
    
    % Using Formula 2: U = lambda / (lambda + mu)
    % data.brlambda is in failures/year
    prob_branch = data.brlambda ./ (data.brlambda + branch_mu);

    %% 3. Combine Results
    % Concatenate generator and branch probabilities into a single vector
    total_unavailability = horzcat(prob_gen, prob_branch);

end
