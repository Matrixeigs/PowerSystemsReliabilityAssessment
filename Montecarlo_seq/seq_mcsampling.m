function stadur = seq_mcsampling(meantime, Ng, Nl, SAMPLE_YEAR, YEAR_HOUR)
%%-------------------------------
%%tic	

	sampleno = SAMPLE_YEAR*YEAR_HOUR;
	stadur = sparse(Ng+Nl, sampleno);%%an Ng+Nl-by-sampleno all zero matrix

	for ieqrow = 1 : Ng+Nl
		
		starow = sparse(1, sampleno);
		staflag = true;
		pointer = 0;
		
		%%-------------------------Monte Carlo state duration sampling---------------%%
		while pointer < sampleno

			if staflag				
				t = floor(-meantime(ieqrow, 1) * log(rand(1)));
				pointer = pointer + t;
			else				
				t = ceil(-meantime(ieqrow, 2) * log(rand(1)));
				starow(1, pointer+1:pointer+t) = 1;
				pointer = pointer + t;
			end

			staflag = ~ staflag;

        end
        
        if size(starow, 2) > sampleno
			starow(:, sampleno+1:end) = [];
		end

        stadur(ieqrow, :) = starow;

		%%---------------------------------------------------------------------------%%
	end

%%	stadur(:, sampleno+1:pointer) = [];
%%	stadur = double(stadur);


%%toc

return
%%----------------------------------end----------------------------------%%












