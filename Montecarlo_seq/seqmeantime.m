% filepath: 
function reliability_data = seqmeantime()
%% ========================================================================
%  COMPONENT RELIABILITY PARAMETERS (MTTF & MTTR)
%  ========================================================================
%  Purpose:
%  Aggregates Mean Time To Failure (MTTF) and Mean Time To Repair (MTTR)
%  for all system components (Generators and Branches).
%
%  Output:
%  - reliability_data: Matrix [(Ng+Nl) x 2]
%       Column 1: MTTF (Hours)
%       Column 2: MTTR (Hours)
%  ========================================================================

    % Load raw reliability data
    raw_data = case24_failrate();

    %% 1. Generator Data
    % Already provided as MTTF and MTTR in hours
    gen_mttf = raw_data.genmttf';
    gen_mttr = raw_data.genmttr';

    %% 2. Branch Data
    % Provided as Failure Rate (lambda, failures/yr) and Repair Duration (r, hours)
    % Convert Lambda to MTTF: MTTF = 8760 / Lambda
    branch_mttf = (8760 ./ raw_data.brlambda)';
    branch_mttr = raw_data.brdur';

    %% 3. Combine into Single Matrix
    % Row indices 1..Ng: Generators
    % Row indices Ng+1..Ng+Nl: Branches
    reliability_data = [
        gen_mttf,    gen_mttr;
        branch_mttf, branch_mttr
    ];

end