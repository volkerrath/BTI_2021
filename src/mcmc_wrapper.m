function []=mcmc_wrapper(job,NAME,DATA,OPTIONS,PARAMS,MODEL)

% SET RANDOM GENERATOR
rng('shuffle');
% randn('state',sum(100*clock));

NAMEI=strcat(['./',NAME,'_',num2str(randi(100000,1,1)),'_job',num2str(job)]);
disp(['>>>>>> ',NAMEI ]);
% % RUN THE CHAIN
filename=strcat(['./',NAMEI,'_save']);OPTIONS.savename=filename;
[results,parchain,reschain,calchain,s2chain] = mcmcrun2(MODEL,DATA,PARAMS,OPTIONS);
filename=strcat(['./',NAMEI,'_final']);save(filename);
% %
end
