%%function [] = seqMain()
% [] = seqMain()
% 
% Monte Carlo sequential simulation main function
% 
% 
% Inputs(void)
% 
% Outputs(void)
% 
% 
% 

%%---------------------------------Commence function-------------------------------------------------------------------------------------------%%

    clc;
    close all;
    clear;
tic

%%	-------------------------------Establish variables and set initial value-------------------------------------------------------------------------------------------%%
	
	%%	testsystem information
	Testsys = loadcase('case24_ieee_rts');
	Nb = size(Testsys.bus, 1);
	Ng = size(Testsys.gen, 1);
	Nl = size(Testsys.branch, 1);

	%%	simulation constant
	YEAR_HOUR = 8736;
	SAMPLE_YEAR = 4000;
	CVEENSTHR = 0.05;

	%%	record error state, hour, load
	errorpointer = 0;
%	errorflag = 0;
	errorstate = zeros(Ng+Nl+3, 100);

	%%	record failure state, hour, load that curtailment is greater than threshold
	curtailpointer = 0;
	CURTAILTHR = 0.01;
	curtailstate = zeros(Ng+Nl+4, 500);


	%%	basic indices for each year
	plc_yr = zeros(SAMPLE_YEAR, 1);
	nlc = zeros(SAMPLE_YEAR, 1);
	dlc = zeros(SAMPLE_YEAR, 1);
	lc = zeros(SAMPLE_YEAR, 1);
	dns = zeros(SAMPLE_YEAR, 1);
	ens = zeros(SAMPLE_YEAR, 1);
    bpii_yr = zeros(SAMPLE_YEAR, 1);

    %%	basic expected indices
	plc = zeros(SAMPLE_YEAR, 1);
	enlc = zeros(SAMPLE_YEAR, 1);
	edlc = zeros(SAMPLE_YEAR, 1);
%	adlc = zeros(SAMPLE_YEAR, 1);
	elc = zeros(SAMPLE_YEAR, 1);
	edns = zeros(SAMPLE_YEAR, 1);
	eens = zeros(SAMPLE_YEAR, 1);
	cveens = zeros(SAMPLE_YEAR, 1);	% coefficient of varience for EENS, one of the criteria for stopping loop
	bpii = zeros(SAMPLE_YEAR, 1);
%	bpeci = zeros(SAMPLE_YEAR, 1);
%	bpaci = zeros(SAMPLE_YEAR, 1);
%	mbpci = zeros(SAMPLE_YEAR, 1);
%	si = zeros(SAMPLE_YEAR, 1);
	
%%	-------------------------------Establish annual load duration curve---------------------------------------------------------------%%
	[busPd, busQd, loadpropor] = anloducurve(YEAR_HOUR); 
	meantime = seqmeantime;

%%	-------------------------MC sequential sampling and identify contingency hour---------------------------------------------------------------%%

%{	
	loadbus = find(Testsys.bus(:, 3) ~= 0);
	sizeloadbus = size(loadbus, 1);
	
%	Testsys.bus(loadbus, 3) = ctbusPd;
%	Testsys.bus(loadbus, 4) = ctbusQd;

	Testsys.Aload = sum(Testsys.bus(:, 3));
	Testsys.gencost = repmat([2 0 0 3 0 0 0 ], Ng, 1);


	Testsys.gencost = vertcat(Testsys.gencost, repmat([2 0 0 3 0 1 0], sizeloadbus, 1));
	
	addgen = Testsys.gen(1:sizeloadbus, :);	

	addgen(:, 1:10) = horzcat(	Testsys.bus(loadbus, 1),					- Testsys.bus(loadbus, 3),	...
								- Testsys.bus(loadbus, 4),				zeros(sizeloadbus, 1),		...
								- Testsys.bus(loadbus, 4),				zeros(sizeloadbus, 1),		...
								Testsys.baseMVA * ones(sizeloadbus, 1),	ones(sizeloadbus, 1),		...	
								zeros(sizeloadbus, 1),					- Testsys.bus(loadbus, 3)	);

	Testsys.gen = vertcat(Testsys.gen, addgen);

	Testsys.bus(loadbus, 3:4) = 0;
%}	
	%	Dispatchable load for objective function of minimum load curtailment
	[Testsys, loadbus]= dispaload(Testsys);

%%	-------------------------------Set Matpower parameters and preliminary test---------------------------------------------------------------%%
	mpopt = mpoption('PF_DC',1,'VERBOSE',0,'OUT_ALL',0,'OPF_ALG_DC',200,'OPF_FLOW_LIM',1);
%    OPFresult = runopf(Testsys, mpopt);

%%	-------------------------------Start loop simulation---------------------------------------------------------------------------%%
	for	iYear = 1 : SAMPLE_YEAR

%%		-------------------------MC sequential sampling and identify contingency hour---------------------------------------------------------------%%       
		stadur = seq_mcsampling(meantime, Ng, Nl, 1, YEAR_HOUR);

		cthour = find(sum(stadur, 1) > 0);
		ctstadur = stadur(:, cthour);
		% ctbusPd = busPd(:, cthour);
		% ctbusQd = busQd(:, cthour);
%%		-------------------------Analyze states in each contingency hour---------------------------------------------------------------%%
		ctana = zeros(2, size(cthour, 2));
        errorflag = zeros(1, size(cthour, 2));
		parfor jHour = 1 : size(cthour, 2)
			try
				ctana(1, jHour) = seq_mcsimulation(ctstadur(:, jHour), loadpropor(1, cthour(1, jHour)) , Testsys, mpopt, loadbus, Ng, Nl);	%% do simulations for each contingency hour
            catch
 %              break;
% 				errorstate(1:Ng+Nl, errorpointer) = ctstadur(:, jHour);
% 				errorstate(Ng+Nl+1, errorpointer) = cthour(1, jHour);
% 				errorstate(Ng+Nl+2, errorpointer) = loadpropor(1, cthour(1, jHour));
% 				errorstate(Ng+Nl+3, errorpointer) = iYear;
    			warning('Problem using Matpower runopf. Assigning a value of 0 to curtailment');
                ctana(1, jHour) = 0;
                errorflag(1, jHour) = 1;
			end
		end		
		ctana(2, :) = double(ctana(1, :) ~= 0);	%% contingency flag
        
%%		---------------------------Record error state that cause Matpower runopf malfunctioned-------------------------------%%      
        numerror = sum(errorflag);
        if numerror
            errorstate(1:Ng+Nl, errorpointer+1:errorpointer+numerror) = ctstadur(:, errorflag == 1 );
            errorstate(Ng+Nl+1, errorpointer+1:errorpointer+numerror) = cthour(1, errorflag == 1);
            errorstate(Ng+Nl+2, errorpointer+1:errorpointer+numerror) = loadpropor(1, cthour(1, errorflag == 1));
            errorstate(Ng+Nl+3, errorpointer+1:errorpointer+numerror) = iYear;
            errorpointer = errorpointer + numerror;
        end
 			
%%		---------------------------Record failure state, hour, load that curtailment is greater than threshold-------------------------------%%
		recordcurtail = find(ctana(1, :) > CURTAILTHR);
		numcurtail = size(recordcurtail, 2);
		if numcurtail > 0
			curtailstate(1:Ng+Nl, curtailpointer+1:curtailpointer+numcurtail) = ctstadur(:, recordcurtail);
			curtailstate(Ng+Nl+1, curtailpointer+1:curtailpointer+numcurtail) = cthour(:, recordcurtail);
			curtailstate(Ng+Nl+2, curtailpointer+1:curtailpointer+numcurtail) = loadpropor(:, recordcurtail);
			curtailstate(Ng+Nl+3, curtailpointer+1:curtailpointer+numcurtail) = iYear;
            curtailstate(Ng+Nl+4, curtailpointer+1:curtailpointer+numcurtail) = ctana(1, recordcurtail);
			curtailpointer = curtailpointer + numcurtail;
        end

%%		-----------------------Calculate Basic Indices per year---------------------------------------------------------------%%	
		plc_yr(iYear, 1) = sum(ctana(2, :) / YEAR_HOUR);
		nlc(iYear, 1) = calnlc(ctana(2, :));		
		dlc(iYear, 1) = sum(ctana(2, :));
		lc(iYear, 1) =	sum(ctana(1, :) * 1);	% load curtailments(MW/yr) LC = sigma CiFi
		dns(iYear, 1) = sum(ctana(1, :) * 1) / YEAR_HOUR; 
		ens(iYear, 1) = sum(ctana(1, :) * 1 * 1);	%   ENS = sigma CiFiDi
		bpii_yr(iYear, 1) = sum(ctana(1, :) * 1 / Testsys.Pload);
		
%%		-----------------------Calculate Basic Expected Indices------------------------------------------------------------------------------------------%%
		plc(iYear, 1) = mean(plc_yr(1:iYear, 1));
		enlc(iYear, 1) = mean(nlc(1:iYear, 1));
 		edlc(iYear, 1) = mean(dlc(1:iYear, 1));
 		elc(iYear, 1) = mean(lc(1:iYear, 1));
 		edns(iYear, 1) = mean(dns(1:iYear, 1));
		
		eens(iYear, 1) = mean(ens(1:iYear, 1));
		cveens(iYear, 1) = std(eens(1:iYear, 1)) / mean(eens(1:iYear));
        plot(eens(1:iYear, 1));
		drawnow;

		bpii = mean(bpii_yr(1:iYear, 1));	

%%		-----------------------Determin if cveens smaller than setting threshold------------------------------------------------------------------------------------------%%		
		if (cveens(iYear, 1) < CVEENSTHR) && (cveens(iYear, 1) ~= 0)
%			disp('Simulation terminated due to achieving accetable accuracy');
%			flag_break = 1;
			break;	% if already achieve acceptable accuracy, break year simulation
		end
        
%%		-----------------------------------------------------------------------------------------------------------------%%
	
	end

%%	-------------------------------Delete all zero placeholders----------------------------------------------------------------------------------%%		
	if iYear < SAMPLE_YEAR
		plc(iYear+1:SAMPLE_YEAR, :) = [];
		enlc(iYear+1:SAMPLE_YEAR, :) = [];
		edlc(iYear+1:SAMPLE_YEAR, :) = [];
		elc(iYear+1:SAMPLE_YEAR, :) = [];
		edns(iYear+1:SAMPLE_YEAR, :) = [];
		eens(iYear+1:SAMPLE_YEAR, :) = [];
		bpii_yr(iYear+1:SAMPLE_YEAR, :) = [];
	end

%%		-----------------------Calculate Remaining Expected Indices------------------------------------------------------------------------------------------%%
	adlc = edlc ./ enlc;
	bpeci = eens / Testsys.Pload;
	bpaci = elc ./ enlc;
	mbpci = edns / Testsys.Pload;
	si = bpeci * 60;
	
toc
	
%%end
%%---------------------------------------------End-------------------------------------------------------------------------%%
