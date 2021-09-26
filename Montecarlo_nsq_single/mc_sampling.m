function eqstatus = mc_sampling(totalprob, SIMUNIT, Ng, Nl)
	
%%	Compare random matrix with augmented failure probability matrix to get the equipment status
%%	1 - failure, 0 - normal

	eqstatus = rand(SIMUNIT, Ng+Nl) < repmat(totalprob, SIMUNIT, 1);
	
	eqstatus(:, 15) = 0;	% neglect the failure probability of synchronous compensator @ bus 14

	eqstatus = sparse(double(eqstatus));

return

%%----------------------------------End-------------------------------------------%%


