function totalprob = failprob()

	failrate = case24_failrate;
%%----------------------------------Generator data----------------------------------%%
%	Compute the failure probability of generators

	probgen = failrate.genmttr ./ (failrate.genmttf + failrate.genmttr);

%{
%	when taking considerations of planned maintainances
	%	componet of maintainance 
	weekhours = 7*24;
	
	mttrp = failrate.genweeks * weekhours;
	genmiup = 8760./ mttrp;
	mttfp = 8760 - mttrp;
	genlambdap = 8760./mttfp;

	%	compoent of forced outage 
	genmiu = 8760./ failrate.genmttr;
	genlambda = 8760./ failrate.genmttf;

%	combination of two parts

	probgen = 	(genlambda * genmiup + genlambdap * genmiu) ./ ...
				(genlambda * genmiup + genlambdap * genmiu + genmiup * genmiu);
%}

%%----------------------------------Branch data----------------------------------%%

	brmiu = 8760 ./ failrate.brdur;
	probbr = failrate.brlambda ./ (failrate.brlambda + brmiu);

%%----------------------------------Combination matrix----------------------------------%%

	totalprob = horzcat(probgen, probbr);

return;
%%----------------------------------End-------------------------------------------%%
