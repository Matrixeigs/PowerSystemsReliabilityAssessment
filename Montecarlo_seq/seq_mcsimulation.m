function curtailment = seq_mcsimulation(ctstadur, loadpropor, Testsys, mpopt, loadbus, Ng, Nl)


%%-------------------------------Load modification--------------------------------------------------------%%

%{	
    genbus = find(Testsys.bus(:, 3) ~= 0);
	sizegenbus = size(genbus, 1);
	
	Testsys.bus(genbus, 3) = ctbusPd;
	Testsys.bus(genbus, 4) = ctbusQd;

	Testsys.Pload = sum(Testsys.bus(:, 3));
	Testsys.gencost = repmat([2 0 0 3 0 0 0 ], Ng, 1);
	

	Testsys.gencost = vertcat(Testsys.gencost, repmat([2 0 0 3 0 1 0], sizegenbus, 1));
	
	addgen = Testsys.gen(1:sizegenbus, :);	

	addgen(:, 1:10) = horzcat(	Testsys.bus(genbus, 1),					- Testsys.bus(genbus, 3),	...
								- Testsys.bus(genbus, 4),				zeros(sizegenbus, 1),		...
								- Testsys.bus(genbus, 4),				zeros(sizegenbus, 1),		...
								Testsys.baseMVA * ones(sizegenbus, 1),	ones(sizegenbus, 1),		...	
								zeros(sizegenbus, 1),					- Testsys.bus(genbus, 3)	);

	Testsys.gen = [Testsys.gen; addgen];

	Testsys.bus(genbus, 3:4) = 0;
	
%}
	sizeloadbus = size(loadbus, 1);
	Testsys.Pload = Testsys.Pload * loadpropor;
	Testsys.gen(Ng+1:Ng+sizeloadbus, [2 3 5 10]) = Testsys.gen(Ng+1:Ng+sizeloadbus, [2 3 5 10]) * loadpropor;

%%-------------------------------Generation modification---------------------------------------------------%%


    Testsys.gen(1:Ng, 8) = ~ ctstadur(1:Ng, 1);
    Testsys.branch(1:Nl, 11) = ~ ctstadur(Ng+1:Ng+Nl, 1);


%%-------------------------------Generation modification------------------------------------------------%%

	Result = runopf(Testsys, mpopt);
	curtailment = Result.f + Testsys.Pload;
	
	if curtailment < 0.1
		curtailment = 0;	
	end

return