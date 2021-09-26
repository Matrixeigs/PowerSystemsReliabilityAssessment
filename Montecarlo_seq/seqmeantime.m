function meantime = seqmeantime()
% meantime = seqmeantime()
% 
% Obtain MTTF and MTTR for each equipment
% 
% Inputs(void)
% 
% Outputs:
% 	meantime	matrix(Ng+Nl, [MTTF, MTTR])	MTTF and MTTR table for each equipment

%{
			equipment	|	MTTF		MTTR	|
			-----------------------------------	|
				1		|			|			|
generator		...		|			|			|		
				Ng		|			|			|
			-----------------------------------	|
				Ng+1	|			|			|
branch			...		|			|			|
				Ng+Nl	|			|			|		
			-----------------------------------	|

%}

%%	---------------------------------Load case file--------------------------------------------------%%
	failrate = case24_failrate;
%	meantime = zeros(Ng+Nl, 2);

%%	----------------------------------Generator data----------------------------------%%	
	genmttf = failrate.genmttf';
	genmttr = failrate.genmttr';

%%	----------------------------------Branch data----------------------------------%%
	brmttf = (8760 ./ failrate.brlambda)';
	brmttr = failrate.brdur';

%%	----------------------------------Aggregate data to form meantime matrix----------------------------------%%

	meantime = [	genmttf,	genmttr;
					brmttf,		brmttr		]; 


return
%%----------------------------------end----------------------------------%%