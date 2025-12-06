function LoadProfile = case24_loadprofile()
%% ========================================================================
%  IEEE RTS-24 LOAD PROFILE DATA
%  ========================================================================
%  Purpose:
%  Provides the hierarchical load factors for the IEEE Reliability Test System.
%  Load(t) = Peak * Weekly(w) * Daily(d) * Hourly(h)
%
%  Data Structure:
%  - MW, MVAr: System Peak Load
%  - weekly:   [1x52] Factors for each week of the year
%  - daily:    [1x7]  Factors for Mon-Sun
%  - hourly:   [24x6] Factors for 24 hours across 3 seasons & 2 day types
%  - busload:  [Nx3]  Bus ID, Active Peak, Reactive Peak
%  ========================================================================

    %% 1. System Peak
    LoadProfile.MW = 2850;
    LoadProfile.MVAr = 580;	

    %% 2. Weekly Factors (52 Weeks)
    % Represents seasonal variation (Winter peak, Summer peak, etc.)
    LoadProfile.weekly = [	
        0.862, 0.900, 0.878, 0.834, ... % Weeks 1-4 (Winter)
        0.880, 0.841, 0.832, 0.806, ... % Weeks 5-8
        0.740, 0.737, 0.715, 0.727, ... % Weeks 9-12 (Spring)
        0.704, 0.750, 0.721, 0.800, ... % Weeks 13-16
        0.754, 0.837, 0.870, 0.880, ... % Weeks 17-20 (Summer start)
        0.856, 0.811, 0.900, 0.887, ... % Weeks 21-24
        0.896, 0.861, 0.755, 0.816, ... % Weeks 25-28
        0.801, 0.880, 0.722, 0.776, ... % Weeks 29-32
        0.800, 0.729, 0.726, 0.705, ... % Weeks 33-36 (Fall)
        0.780, 0.695, 0.724, 0.723, ... % Weeks 37-40
        0.743, 0.744, 0.800, 0.881, ... % Weeks 41-44
        0.885, 0.909, 0.940, 0.890, ... % Weeks 45-48 (Winter)
        0.942, 0.970, 1.000, 0.952      % Weeks 49-52 (Peak at Wk 51)
    ];

    %% 3. Daily Factors (Mon-Sun)
    %                       Mon   Tue   Wed   Thu   Fri   Sat   Sun	
    LoadProfile.daily = [	0.93, 1.00, 0.98, 0.96, 0.94, 0.77, 0.75 ];

    %% 4. Hourly Factors (24 Hours)
    % Columns:
    % 1: Winter Weekday, 2: Winter Weekend
    % 3: Summer Weekday, 4: Summer Weekend
    % 5: Spr/Fall Wkdy,  6: Spr/Fall Wknd
    LoadProfile.hourly = [	
        0.67, 0.78, 0.64, 0.74, 0.63, 0.75; % Hour 1 (12am-1am)
        0.63, 0.72, 0.60, 0.70, 0.62, 0.73; % Hour 2
        0.60, 0.68, 0.58, 0.66, 0.60, 0.69; % Hour 3
        0.59, 0.66, 0.56, 0.65, 0.58, 0.66; % Hour 4
        0.59, 0.64, 0.56, 0.64, 0.59, 0.65; % Hour 5
        0.60, 0.65, 0.58, 0.62, 0.65, 0.65; % Hour 6
        0.74, 0.66, 0.64, 0.62, 0.72, 0.68; % Hour 7
        0.86, 0.70, 0.76, 0.66, 0.85, 0.74; % Hour 8
        0.95, 0.80, 0.87, 0.81, 0.95, 0.83; % Hour 9
        0.96, 0.88, 0.95, 0.86, 0.99, 0.89; % Hour 10
        0.96, 0.90, 0.99, 0.91, 1.00, 0.92; % Hour 11
        0.95, 0.91, 1.00, 0.93, 0.99, 0.94; % Hour 12
        0.95, 0.90, 0.99, 0.93, 0.93, 0.91; % Hour 13
        0.95, 0.88, 1.00, 0.92, 0.92, 0.90; % Hour 14
        0.93, 0.87, 1.00, 0.91, 0.90, 0.90; % Hour 15
        0.94, 0.87, 0.97, 0.91, 0.88, 0.86; % Hour 16
        0.99, 0.91, 0.96, 0.92, 0.90, 0.85; % Hour 17
        1.00, 1.00, 0.96, 0.94, 0.92, 0.88; % Hour 18 (Peak)
        1.00, 0.99, 0.93, 0.95, 0.96, 0.92; % Hour 19
        0.96, 0.97, 0.92, 0.95, 0.98, 1.00; % Hour 20
        0.91, 0.94, 0.92, 1.00, 0.96, 0.97; % Hour 21
        0.83, 0.92, 0.93, 0.93, 0.90, 0.95; % Hour 22
        0.73, 0.87, 0.87, 0.88, 0.80, 0.90; % Hour 23
        0.63, 0.81, 0.72, 0.80, 0.70, 0.85  % Hour 24
    ];

    %% 5. Bus Load Data
    % [ BusID, Active_Peak(MW), Reactive_Peak(MVAr) ]
    LoadProfile.busload = [
        1,  108, 22;	
        2,  97,  20;	
        3,  180, 37;	
        4,  74,  15;	
        5,  71,  14;	
        6,  136, 28;
        7,  125, 25;	
        8,  171, 35;	
        9,  175, 36;	
        10, 195, 40;				
        13, 265, 54;	
        14, 194, 39;	
        15, 317, 64;	
        16, 100, 20;					
        18, 333, 68;	
        19, 181, 37;	
        20, 128, 26
    ];

end

