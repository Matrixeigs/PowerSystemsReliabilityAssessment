function [busPd, busQd, loadpropor] = anloducurve(YEAR_HOUR)
%[busPd, busQd, loadpropor] = anloducurve(YEAR_HOUR)
%
%Establish the hourly load duration curve based on IEEE Reliability Test System 1979.
%Normally, it consists of 8736 hours. Each bus has same load ratios at each certain time.
%
%Inputs:
%	YEAR_HOUR:	scalar						the number of points of annual load duration curve, 
%											it actually should be 8760 (365*24), but IEEE RTS-79 
%											only provides 8736 hours (52*7*24).
%
%Outputs:
%	busPd:		Matrix(Ng+Nl, YEAR_HOUR)	Active power load at each bus for each hour
%	busQd:		Matrix(Ng+Nl, YEAR_HOUR)	Reactive power load at a each bus for each hour
%	loadpropor	Vector(1, YEAR_HOUR)		Load ratio at each bus for each hour

%%%-------------------------------Function start--------------------------------------------------------------------------%%
%tic
	%%-------------------------------Load basic porfile-------------------------%%

	Loadprofile = case24_loadprofile;	
	loadpropor = zeros(1, YEAR_HOUR);

	%%-------------------------------Circulation structure for hourly load curve-------------%%
	for ihour = 1 : YEAR_HOUR
		
		%%--------------------Identify a certain hour to corresponding week------------------------%%	
		%% week classification
		corweek = ceil(ihour/168);
        
		if corweek <= 8 || corweek >=44
			weekflag = 'winter';
		elseif corweek >= 18 && corweek <= 30
			weekflag = 'summer';
		elseif corweek <= 17 || corweek >= 31
			weekflag = 'springfall';
		else
			disp('error');
		end
		%%-----------------------Identify a certain hour to corresponding day------------------------------------------%%			
		
		corday = ceil(mod(ihour/24, 7));
        
        if corday == 0	%%sunday the remainder is 0, needs to be calibrated
            corday = 7; 
        end
		
		%% day classification
		if corday <= 5
			dayflag = 'weekday';
		else
			dayflag = 'weekend';
		end
		%%---------------------Identify a certain hour to corresponding hour-----------------------%%
		corhour = mod(ihour, 24);
        
        if corhour == 0	%%for 23:00-24:00 time interval, the remainder is 0, needs to be calibrated
            corhour = 24;
        end

		%%-----------------------Obtain weekly ratio, daily ratio and hourly ratio ---------------------%%

		weekpropor = Loadprofile.weekly(corweek);
		daypropor = Loadprofile.daily(corday);

		switch weekflag
			case 'winter'
				switch dayflag
					case 'weekday'
						hourpropor = Loadprofile.hourly(corhour, 1);
					case 'weekend'
						hourpropor = Loadprofile.hourly(corhour, 2);
				end
			case 'summer'
				switch dayflag
					case 'weekday'
						hourpropor = Loadprofile.hourly(corhour, 3);
					case 'weekend'
						hourpropor = Loadprofile.hourly(corhour, 4);
				end
			case 'springfall'
				switch dayflag
					case 'weekday'
						hourpropor = Loadprofile.hourly(corhour, 5);
					case 'weekend'
						hourpropor = Loadprofile.hourly(corhour, 6);
				end
		end

		%%-----------------------Obtain hourly load duration curve ---------------------%%
		loadpropor(1, ihour) = weekpropor * daypropor * hourpropor;			
    end

	%%-------------------------------Obtain load duration curve for each bus---------------------%%    
	
	busPd = Loadprofile.busload(:, 2) * loadpropor;
	busQd = Loadprofile.busload(:, 3) * loadpropor;
    
%toc   
return
%%-------------------------------function end---------------------%%		
