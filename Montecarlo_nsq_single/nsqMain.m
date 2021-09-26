%%function [] = nsqMain()
%{

0


System Indices:

	ENLC	Expected Number of Load Curtailments (occ./yr)
	EDLC	Expected Duration of Load Curtailments (hr/yr)
	PLC		Probability of Load Curtailments
	EDNS  	Expected Demand Not Supplied (MW)
	EENS 	Expected Energy Not Supplied (MWh/yr)
	BPII	Bulk Power Interruption Index (MW/MW-yr)
	BPECI	Bulk Power/Energy Curtailment Index (MWh/MW-yr)
	BPACI	Bulk Power-supply Average MW Curtailment Index (MW/disturbance)
	MBPCI	Modified Bulk Power Curtailment Index (MW/MW)
	SI 		Severity Index (system minutes/yr)     



%}

	clear
	close all
	clc

tic

%%	Load test system
	
	Testsys = loadcase('case24_ieee_rts');
	Nb = size(Testsys.bus, 1);
	Ng = size(Testsys.gen, 1);
	Nl = size(Testsys.branch, 1);


%%	Set initial value

	BETAlimit = 0.0017;
	ITER_max = 100000;
	SIMUNIT = 100;
	

	iter = 0;
	lole = 0;
	edns = 0;
	betavalue = inf;
	row_index = 0;

%% Build matrices that have fix dimension to avoid changing size in each loop;
	eqstatus_total = sparse(ITER_max, Ng+Nl+3);
	beta_table = zeros(1, ITER_max / SIMUNIT);
	edns_table = zeros(1, ITER_max / SIMUNIT);
	lole_table = zeros(1, ITER_max / SIMUNIT);
    plc_table = zeros(1, ITER_max / SIMUNIT);

% 	eqstatus_total = [];
% 	beta_table = [];
% 	eens_table = [];
% 	lolp_table = [];


%%	Dispachable load model, enable different control method
%%	Change runopf objective function to minumum load curtailment through changing Pg , Ld, Gcost;
%%	
%{
t1 = tic
	Testsys.load = sum(Testsys.bus(:, 3));
	Testsys.gencost = repmat([2 0 0 3 0 0 0], Ng, 1);

	for i=1 : Nb
	    if Testsys.bus(i, 3) ~= 0
	        Testsys.gencost =[Testsys.gencost; [2,0,0,3,0,1,0]];%%负荷等级更改！
	        index = zeros(1, 10);
	        index(1) = Testsys.bus(i, 1);	% Bus Index
	        index(2) = -Testsys.bus(i, 3);	% Active Power
	        index(3) = -Testsys.bus(i, 4);	% Reactive Power
	        index(4) = 0;
	        index(5) = -Testsys.bus(i,4);	% Minimum Reactive Power
	        index(7) = Testsys.baseMVA;
	        index(8) = 1;	% On Line Or Off Line
	        index(9) = 0;
	        index(10) = -Testsys.bus(i,3);	% Minimum Active Power
	        Index = Testsys.gen(1, :);
	        Index(1: 10) = index;
	        Testsys.gen = [Testsys.gen;Index];
	        clear Index index;
	        %Remove Load
	        Testsys.bus(i,3)=0;
	        Testsys.bus(i,4)=0;
	    end
	end
%}
%	Add controllable load procedure
%%{
%t1 = tic	
	genbus = find(Testsys.bus(:, 3) ~= 0);
	sizegenbus = size(genbus, 1);

	Testsys.load = sum(Testsys.bus(:, 3));
	Testsys.gencost = repmat([2 0 0 3 0 0 0 ], Ng, 1);%%all gen cost is 0

	
	
    %% treat all load as negtive generator and set their parameters, then add these vitual generators to real gens
	Testsys.gencost = vertcat(Testsys.gencost, repmat([2 0 0 3 0 1 0], sizegenbus, 1));%%treat load as gen
	
	Index = Testsys.gen(1:sizegenbus, :);	

	Index(:, 1:10) = horzcat(	Testsys.bus(genbus, 1),					- Testsys.bus(genbus, 3),	...
								- Testsys.bus(genbus, 4),				zeros(sizegenbus, 1),		...
								- Testsys.bus(genbus, 4),				zeros(sizegenbus, 1),		...
								Testsys.baseMVA * ones(sizegenbus, 1),	ones(sizegenbus, 1),		...	
								zeros(sizegenbus, 1),					- Testsys.bus(genbus, 3)	);

	Testsys.gen = [Testsys.gen; Index];

	clear Index;

	Testsys.bus(genbus, 3:4) = 0;%% set original load to 0 at their buses
%%}	

%toc(t1)

%%	Load equipment failure probability
	totalprob = failprob; %%failure probability(gen,branch)


	mpopt = mpoption('PF_DC',1,'VERBOSE',0,'OUT_ALL',0,'OPF_ALG_DC',200,'OPF_FLOW_LIM',1);

	result = runopf(Testsys, mpopt);    

	while betavalue > BETAlimit && iter < ITER_max;
%tstart1 = tic;
	%%----------------------------------Monte Carlo sampling---------------------%%	

		eqstatus_indi = mc_sampling(totalprob, SIMUNIT, Ng, Nl);
        
        eqstatus_indi(:, Ng+Nl+1) = 1;
        
        eqstatus_indi(:, Ng+Nl+2:Ng+Nl+3) = 0;

		[eqstatus_indi, ia1, ic1] = unique(eqstatus_indi, 'rows', 'stable'); % check the repeatd sample in eqstatus_indi;

		for i = 1 : size(eqstatus_indi, 1)
			eqstatus_indi(i, Ng+Nl+1) = sum(double(ic1 == i));
		end
%tend1 = toc(tstart1)
	%%	Obtain a sampling status matrix that already combines all the repeated rows and set the account results to the column Ng+Nl+1;

	%%----------------------------------Monte Carlo simulation--------------------%%
	%%	Obtain OPF results for each row and set the results to the column NG+Nl+2, set the load shedding sign to the column Ng+Nl+3;

		if iter 
	%%	Combine the new eqstatus_indi to eqstatus_total, repeated rows need to be account to the column Ng+Nl+1;
	%%	Then compute the OPF for new status in eqstatus_indi; 				
%tstart2 = tic;
			
			[eqstatus_temp, ia2, ib2] = intersect(eqstatus_total(1:row_index, 1:Ng+Nl), eqstatus_indi(:, 1:Ng+Nl), 'rows', 'stable'); % check the repeated status 

			eqstatus_total(ia2, Ng+Nl+1) = eqstatus_total(ia2, Ng+Nl+1) + eqstatus_indi(ib2, Ng+Nl+1); % accumulating count

			eqstatus_indi(ib2, :) = []; % Delete repeted rows in eqstatus_indi;
			
			
			%{
			for irow = 1 : size(eqstatus_indi, 1)
				for jcolumn = 1 : size(eqstatus_total, 1)
					if isequal(eqstatus_indi(irow, 1:Ng+Nl), eqstatus_total(jcolumn, 1:Ng+Nl))
						eqstatus_total(jcolumn, Ng+Nl+1) = eqstatus_total(jcolumn, Ng+Nl+1) + eqstatus_indi(irow, Ng+Nl+1);
						eqstatus_indi(irow, Ng+Nl+1) = 0;
						break
					end
				end
			end

			zerorow = find(eqstatus_indi(:, Ng+Nl+1) == 0);
			eqstatus_indi(zerorow, :) = [];

			clear zerorow;

			%}


%tend2 = toc(tstart2)            
            parfortemp = zeros(size(eqstatus_indi, 1), 2);

%tstart3 = tic;
			parfor i = 1:size(eqstatus_indi, 1)
				
                parfortemp(i, 1) = mc_simulation(eqstatus_indi(i, 1:Ng+Nl), Testsys, mpopt, Ng, Nl);
				
            end
 %tend3 = toc(tstart3)  
 %tstart4 = tic;       
            parfortemp(:, 2) = double(parfortemp(:, 1) ~= 0);
            
            eqstatus_indi(:, Ng+Nl+2:Ng+Nl+3) = parfortemp;

%            eqstatus_total = [eqstatus_total; eqstatus_indi];

			eqstatus_total(row_index+1:row_index+size(eqstatus_indi), :) = eqstatus_indi;

			row_index = row_index + size(eqstatus_indi, 1); % Accumulate row_index
%tend4 = toc(tstart4)
        else
            
            parfortemp = zeros(size(eqstatus_indi, 1), 2);
            % parfortemp first column is how much load curtailed in state i
            % parfortemp second column is whether load is curtailed
            
			parfor i = 1:size(eqstatus_indi, 1)
				parfortemp(i, 1) = mc_simulation(eqstatus_indi(i, 1:Ng+Nl), Testsys, mpopt, Ng, Nl);
            end
				
			
            parfortemp(:, 2) = double(parfortemp(:, 1) ~= 0);
            
            eqstatus_indi(:, Ng+Nl+2:Ng+Nl+3) = parfortemp;

%            eqstatus_total = eqstatus_indi;

			eqstatus_total(row_index+1:row_index+size(eqstatus_indi), :) = eqstatus_indi;

			row_index = row_index + size(eqstatus_indi, 1); % Accumulate row_index 



		end

	%%----------------------------------Indices computation--------------------------%%
%tstart5 = tic;
		edns = 	sum( eqstatus_total(1:row_index, Ng+Nl+1) .* eqstatus_total(1:row_index, Ng+Nl+2))	...
				/ (iter + SIMUNIT) ;
		lole = 	sum( eqstatus_total(1:row_index, Ng+Nl+1) .* eqstatus_total(1:row_index, Ng+Nl+3))	...
				/ (iter + SIMUNIT) * 8760 ;
         plc= sum( eqstatus_total(1:row_index, Ng+Nl+1) .* eqstatus_total(1:row_index, Ng+Nl+3))	...
            / (iter + SIMUNIT);
		betavalue = sqrt(sum(eqstatus_total(1:row_index, Ng+Nl+1) .* (eqstatus_total(1:row_index, Ng+Nl+2) - edns) .^2))	...
					/ (iter + SIMUNIT) / edns;


% 		eens = 	sum( eqstatus_total(:, Ng+Nl+1) .* eqstatus_total(:, Ng+Nl+2))	...
% 				/ (iter + SIMUNIT) ;
% 		lolp = 	sum( eqstatus_total(:, Ng+Nl+1) .* eqstatus_total(:, Ng+Nl+3))	...
% 				/ (iter + SIMUNIT) ;
% 		betavalue = sqrt(sum(eqstatus_total(:, Ng+Nl+1) .* (eqstatus_total(:, Ng+Nl+2) - eens) .^2))	...
% 					/ (iter + SIMUNIT) / eens;

	%%	Record beta, eens, lolp values of each iteration;
		beta_table(1, (iter + SIMUNIT)/SIMUNIT) = betavalue;
		edns_table(1, (iter + SIMUNIT)/SIMUNIT) = edns;
		lole_table(1, (iter + SIMUNIT)/SIMUNIT) = lole;
         plc_table(1, (iter + SIMUNIT)/SIMUNIT) = plc;


% 		beta_table = [beta_table, betavalue];
% 		eens_table = [eens_table, eens];
% 		lolp_table = [lolp_table, lolp];

	%%	incremental value to iter;
	
		iter = iter + SIMUNIT;
	
	%% 	plot the beta value		

	%%	plot(1:iter/SIMUNIT, beta_table(1:iter/SIMUNIT)); 
	%%	drawnow;
%tend5 = toc(tstart5)
	end

toc

%%end
%%----------------------------------nsqMain end------------------------------%%


