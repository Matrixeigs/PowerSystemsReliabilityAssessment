function [Testsys, loadbus] = dispaload(Testsys)
	
    Nb = size(Testsys.bus, 1);
	Ng = size(Testsys.gen, 1);
	Nl = size(Testsys.branch, 1);

	loadbus = find(Testsys.bus(:, 3) ~= 0);
	sizeloadbus = size(loadbus, 1);

	Testsys.Pload = sum(Testsys.bus(:, 3));
	Testsys.gencost = repmat([2 0 0 3 0 0 0 ], Ng, 1);
	
	Testsys.gencost = vertcat(Testsys.gencost, repmat([2 0 0 3 0 1 0], sizeloadbus, 1));
	
% 	addgen = Testsys.gen(1:sizeloadbus, :);	
    addgen = zeros(sizeloadbus, size(Testsys.gen, 2));

	addgen(:, 1:10) = horzcat(	Testsys.bus(loadbus, 1),					- Testsys.bus(loadbus, 3),	...
								- Testsys.bus(loadbus, 4),				zeros(sizeloadbus, 1),		...
								- Testsys.bus(loadbus, 4),				zeros(sizeloadbus, 1),		...
								Testsys.baseMVA * ones(sizeloadbus, 1),	ones(sizeloadbus, 1),		...	
								zeros(sizeloadbus, 1),					- Testsys.bus(loadbus, 3)	);

	Testsys.gen = vertcat(Testsys.gen, addgen);

	Testsys.bus(loadbus, 3:4) = 0;

end
