function []=Fwd_gsth(name)
% GSTH_Tikh performs paleoclimatic inversion in parameter
% space for a single SITE

% important index arrays for indirekt addresses:
% id        points to nodes defining obs
% ip        associates parameter values of diffusivity and production to cells
% it        accociates paleotemperatures to temporal grid cells
%
% structures used:
% sitemod     information defining parameter for given site
% fwdpar      numerical control for forward modeling
% invpar      numerical control for inversion
%
% V. R. March 2019

inv_debug=1;
year2sec=31557600;

load('common');

% PATHS
addpath([srcpath,'/src']);
addpath([srcpath,'/tools']);

dfmt=1;
%ffmt='.zip';
%archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


disp([' ']);
disp(strcat([ '============================================']));
disp(['Nonlinear Thermal Modelling, stationary case'])
disp([name,' --- ',datestr(now,dfmt)])
disp(strcat([ '============================================']));
if run_parallel
    disp(['Number of processes: ',num2str(parcors)])
    disp(strcat([ '============================================']));
end
disp([' ']);


% if run_parallel
%     if isempty(gcp('nocreate')), parpool(parcors); end
% end

% if exist('matlabpool')==0
% if isempty(gcp('nocreate')), parpool(parcors); end
% else
%      if matlabpool('size') == 0
%         matlabpool (parcors)
%      end
% end

% READ PARAMETERS AND DATA FOR SITE
disp([' ']); disp([' ...read & organize obs ' ]);
nobs=0;

F=strcat([name '_SiteMod.mat']);
if exist(F,'file')
    load(F);
    disp([' >>>>> site obs read from: ' F]);
    disp([' ']);
    mstruct(sitemod);
else
    error([F ' does not exist in path!']);
end


% SITE SPECIFIC PATHS
addpath([srcpath,'/src/props/',props]);
disp([' ']); disp([' ...local property model: ' props]);

% READ LOCAL MODELLING PARAMETERS
F=strcat([name,'_FwdPar.mat']);
if exist(F,'file')
    disp([' ']);
    disp([mfilename,': local inversion pars loaded from file ',F]);
    load(F);mstruct(fwdpar);
else
    error([F ' does not exist in path! STOP.']);
end


% ITERATION
disp(['   ']);disp([' ... start iteration   ']);

% SOLVE FORWARD PROBLEM

dt=diff(t);nt=length(t);
dz=diff(z);nz=length(z);

    T0=heat1dns(k, kA, kB,h,r,p,qb,GST(1),dz,ip,maxitnl,tolnl,freeze,out);


[Tcalc,dT,Q,K]=heat1dns(k, kA, kB,h,r,p,qb,GST(1),dz,ip,maxitnl,tolnl,freeze,out);


%==========================================================================
%  Postpocessing
%==========================================================================
% if exist('matlabpool')==0
% if ~isempty(gcp('nocreate')),delete(gcp); end
% else
%      if matlabpool('size') ~= 0
%         matlabpool ('close')
%      end
% end

% SAVE DATA
filename=strcat([name,'_FwdStat.mat']);
disp(['   ']);disp([' ...save results to: ',filename]);
save(filename,'Tcalc','z',)

S=fopen('INFO.dat','a+');
sline=strcat([...
    ' ',num2str(rms,'% 10.6g'),...
    ' ',num2str(mae,'%10.6g'),...
    ' ' name]);

fprintf(S,'%s \n',sline);
fclose(S);


clear



clear


disp(['   '])

