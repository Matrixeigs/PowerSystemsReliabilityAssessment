function [busPd, busQd, load_factors] = anloducurve(TOTAL_HOURS)
%% ========================================================================
%  ANNUAL LOAD DURATION CURVE GENERATOR
%  ========================================================================
%  Purpose:
%  Constructs the hourly load profile based on the IEEE RTS-79 specification.
%  The load model is hierarchical:
%  Load(t) = Peak_Load * Weekly_Factor * Daily_Factor * Hourly_Factor
%
%  Inputs:
%  - TOTAL_HOURS: Number of hours to generate (usually 8736 or 8760)
%
%  Outputs:
%  - busPd: Matrix of Active Load [NumBuses x Hours]
%  - busQd: Matrix of Reactive Load [NumBuses x Hours]
%  - load_factors: Vector of scaling factors [1 x Hours]
%  ========================================================================

    % Load raw profile data (percentages)
    ProfileData = case24_loadprofile();
    load_factors = zeros(1, TOTAL_HOURS);

    %% Generate Hourly Factors
    for hour_idx = 1 : TOTAL_HOURS
        
        % 1. Determine Week of Year (1-52)
        week_idx = ceil(hour_idx / 168); % 168 hours per week
        
        % Determine Season based on Week
        if week_idx <= 8 || week_idx >= 44
            season = 'winter';
        elseif week_idx >= 18 && week_idx <= 30
            season = 'summer';
        else
            season = 'springfall';
        end
        
        % 2. Determine Day of Week (1=Mon ... 7=Sun)
        day_idx = ceil(mod(hour_idx/24, 7));
        if day_idx == 0; day_idx = 7; end
        
        if day_idx <= 5
            day_type = 'weekday';
        else
            day_type = 'weekend';
        end
        
        % 3. Determine Hour of Day (1-24)
        hour_of_day = mod(hour_idx, 24);
        if hour_of_day == 0; hour_of_day = 24; end

        % 4. Retrieve Factors
        week_factor = ProfileData.weekly(week_idx);
        day_factor  = ProfileData.daily(day_idx);
        
        % Retrieve Hourly Factor based on Season and Day Type
        % ProfileData.hourly columns:
        % 1: Winter Wkdy, 2: Winter Wknd
        % 3: Summer Wkdy, 4: Summer Wknd
        % 5: Spr/Fall Wkdy, 6: Spr/Fall Wknd
        
        col_idx = 0;
        switch season
            case 'winter'
                if strcmp(day_type, 'weekday')
                    col_idx = 1;
                else
                    col_idx = 2;
                end
            case 'summer'
                if strcmp(day_type, 'weekday')
                    col_idx = 3;
                else
                    col_idx = 4;
                end
            case 'springfall'
                if strcmp(day_type, 'weekday')
                    col_idx = 5;
                else
                    col_idx = 6;
                end
        end
        
        hour_factor = ProfileData.hourly(hour_of_day, col_idx);

        % 5. Calculate Combined Factor
        load_factors(hour_idx) = week_factor * day_factor * hour_factor;
    end

    %% Apply to Bus Loads
    % busload(:,2) is Peak Active Load, busload(:,3) is Peak Reactive Load
    busPd = ProfileData.busload(:, 2) * load_factors;
    busQd = ProfileData.busload(:, 3) * load_factors;

end
