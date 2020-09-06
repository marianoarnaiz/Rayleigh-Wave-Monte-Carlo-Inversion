% function of calculating rayleigh wave dispersion 
% 	model: [thickness vp vs]
% 	vec_T: period 
%   VpVs: Vp/Vs Relationship
% 	data: [phv];
function phv = dispR_surf96(vec_T,model,VpVs)

system('rm start.mod');
system('rm temp.dsp');

% make surf96 par file
make_par_surf96('J');

% write model to file

datatemp(:,1) = vec_T(:);
datatemp(:,2) = 3;
datatemp(:,3)  =0.1;

writedisp_surf96(datatemp,'rc','R','C'); % write temp data into dispersion file 
writedisp_surf96(datatemp,'ru','R','U'); % write temp data into dispersion file 

system('cat rc ru > disp_obs.dsp'); 

writemod_surf96(model,'start.mod',VpVs);

system('surf96 1 27 temp.dsp > NUL');
% read dispersion from tmep file
data=readdisp_surf96('temp.dsp');

if(~isempty(data))
	phv=data(:,2);
else
	phv=[]
end
%
system('surf96 39'); % clean up
system('rm ru rc NUL');
	return
end


