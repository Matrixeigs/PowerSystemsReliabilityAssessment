% filepath: /Users/matrixeigs/Codes/PowerSystemsReliabilityAssessment/Montecarlo_nsq_single/case24_failrate.m
function failrate = case24_failrate()
%% ========================================================================
%  IEEE RTS-24 RELIABILITY DATA
%  ========================================================================
%  Purpose:
%  Provides the reliability parameters (MTTF, MTTR, etc.) for the 
%  IEEE Reliability Test System (RTS-24).
%
%  Data Structure:
%  - genmttf: Generator Mean Time To Failure (Hours)
%  - genmttr: Generator Mean Time To Repair (Hours)
%  - brlambda: Branch Failure Rate (Failures/Year)
%  - brdur:    Branch Average Repair Duration (Hours)
%  ========================================================================

    %% --------------------------------------------------------------------
    %  1. GENERATION DATA
    %  --------------------------------------------------------------------
    
    % Generator MTTF (Mean Time To Failure) in Hours
    % Represents the average time a unit operates before failing
    failrate.genmttf = [
        450,   450,   1960,  1960,  450,   ... % Bus 1, 2
        450,   1960,  1960,  1200,  1200,  ... % Bus 2, 7
        1200,  950,   950,   950,   10000, ... % Bus 7, 13, 14 (Sync Comp)
        2940,  2940,  2940,  2940,  2940,  ... % Bus 15, 16
        960,   960,   1100,  1100,  1980,  ... % Bus 18, 21
        1980,  1980,  1980,  1980,  1980,  ... % Bus 22, 23
        960,   960,   1150                     % Bus 23
    ];

    % Generator MTTR (Mean Time To Repair) in Hours
    % Represents the average time to fix a unit after failure
    failrate.genmttr = [
        50,    50,    40,    40,    50,    ...
        50,    40,    40,    50,    50,    ...
        50,    50,    50,    50,    0.1,   ... % Sync Comp (Bus 14) repairs quickly
        60,    60,    60,    60,    60,    ...
        40,    40,    150,   150,   20,    ...
        20,    20,    20,    20,    20,    ...
        40,    40,    100
    ];

    % Generator Scheduled Maintenance (Weeks/Year)
    % (Currently unused in the basic NSQ model but kept for reference)
    failrate.genweeks = [
        2, 2, 3, 3, 2, ...
        2, 3, 3, 3, 3, ...
        3, 4, 4, 4, 0.1, ...
        2, 2, 2, 2, 2, ...
        4, 4, 6, 6, 2, ...
        2, 2, 2, 2, 2, ...
        4, 4, 5
    ];

    %% --------------------------------------------------------------------
    %  2. TRANSMISSION BRANCH DATA
    %  --------------------------------------------------------------------
    
    % Branch Failure Rate (Lambda) in Failures/Year
    failrate.brlambda = [
        0.24, 0.51, 0.33, 0.39, 0.48, 0.38, ...
        0.02, 0.36, 0.34, 0.33, 0.30, 0.44, ...
        0.44, 0.02, 0.02, 0.02, 0.02, 0.40, ...
        0.39, 0.40, 0.52, 0.49, 0.38, 0.33, ...
        0.41, 0.41, 0.41, 0.35, 0.34, 0.32, ...
        0.54, 0.35, 0.35, 0.38, 0.38, 0.34, ...
        0.34, 0.45
    ];
    
    % Branch Repair Duration (r) in Hours
    failrate.brdur = [
        16, 10, 10, 10, 10, 768, 10, 10, 35, 10, 10, 10, ...
        10, 768, 768, 768, 768, 11, 11, 11, 11, 11, 11, 11, ...
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, ...
        11, 11
    ];

end
