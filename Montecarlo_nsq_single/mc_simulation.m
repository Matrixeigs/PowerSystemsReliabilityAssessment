function dns = mc_simulation(eqstatus, Testsys, mpopt, Ng, Nl)

%%----------------------------------Obtain EENS for each sample----------------------------------------------%%
	
	statusgen = eqstatus(:, 1:Ng);
	Testsys.gen(1:Ng, 8) = ~ statusgen';

	statusbranch = eqstatus(:, Ng+1:Ng+Nl);
	Testsys.branch(1: Nl, 11) = ~ statusbranch';
%%-------------------------------Run OPF calculation------------------------------------------------------%%
	Result = runopf(Testsys, mpopt);
	dns = Result.f + Testsys.load;

	if dns < 0.1
		dns = 0;
	end
	
return
%%----------------------------------End----------------------------------------------%%


