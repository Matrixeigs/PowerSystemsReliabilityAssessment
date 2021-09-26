function nlc = calnlc(sysstatus)
%{
	
 Calculates the number of load curtailment for one year.


 Inputs:
 	sysstatus: the hourly system status.

 Outputs:
 	nlc: number of load curtailments.


For instance, sysstatus = [1 1 1 0 0 0 1 1 0 0 0 0 0 0 1] 
denotes 3 power outages in this year.

0 - system in normal operation.
1 - system in corresponding hour has a blackout.

Power outages refers to each continuous blackouts.



%}
%%-----------------------Calculate the number of load curtailments--------------------------------------------------%%
		
		outagesign = diff(sysstatus);	%% [sysstatus(column+1) - sysstatus(column)]
		transindex = find(outagesign ~= 0);	%% find the position of nonzero numbers

		%%-----------------------------------Get correction term-------------------------------------------------%%

		if ~isempty(transindex)
			%% if the first nonzero element is -1 and the final nonzero element is 1, there need to be corrected to obtain the number of power outages
			if (outagesign(1, transindex(1))) == -1 && (outagesign(1, transindex(end)) == 1) 
				corterm = 1;
			else
				corterm = 0;
			end
		%% otherwise, it doesn't need to be corrected	
		else
			
			corterm = 0;

		end	
		%%-----------------------Get nlc--------------------------------------------------%%
		nlc = ceil(sum(abs(outagesign))/2) + corterm;

return
%%-----------------------calnlc end--------------------------------------------------%%
		
