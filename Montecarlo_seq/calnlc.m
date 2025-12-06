function num_curtailments = calnlc(system_status_series)
%% ========================================================================
%  CALCULATE NUMBER OF LOAD CURTAILMENTS (FREQUENCY)
%  ========================================================================
%  Purpose:
%  Counts the number of distinct failure events in a time series.
%  A "failure event" is defined as a transition from Normal (0) to Failed (1).
%  Continuous hours of failure count as a single event.
%
%  Input:
%  - system_status_series: Binary vector [1 x Hours] (1=Failure, 0=Normal)
%
%  Output:
%  - num_curtailments: Integer count of distinct events
%  ========================================================================

    % Calculate transitions using diff()
    % 1  -> 0 : -1 (Repair)
    % 0  -> 1 : +1 (Failure)
    % 1  -> 1 :  0 (Ongoing Failure)
    % 0  -> 0 :  0 (Ongoing Normal)
    transitions = diff(system_status_series);
    
    % Count the number of +1 transitions (Normal to Failure)
    % This captures the start of every outage event
    failures_starts = sum(transitions == 1);
    
    % Edge Case: If the year STARTS in a failure state (Hour 1 = 1)
    % diff() won't catch this because there's no previous hour.
    if system_status_series(1) == 1
        failures_starts = failures_starts + 1;
    end
    
    num_curtailments = failures_starts;

end

