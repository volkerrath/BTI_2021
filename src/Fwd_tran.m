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
%
% V. R. Oct 2019

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

% READ PARAMETERS AND DATA FOR SITE
disp([' ']); disp([' ...read & organize obs ' ]);
nobs=0;

F=strcat([name '_SiteMod.mat']);
if exist(F,'file')
    load(F);
    disp([' >>>>> site model read from: ' F]);
    disp([' ']);
    mstruct(sitemod);
else
    error([F ' does not exist in path!']);
end

F=strcat([name '_SiteObs.mat']);
if exist(F,'file')
    load(F);
    disp([' >>>>> site obs read from: ' F]);
    disp([' ']);
    mstruct(siteobs);
else
    Tobs=[];
    disp([' >>>>> no site obs found! ']);
    disp([' ']);
end

% READ LOCAL MODELLING PARAMETERS
F=strcat([name,'_FwdPar.mat']);
if exist(F,'file')
    disp([' ']);
    disp([mfilename,': local modelling pars loaded from file ',F]);
    load(F);
    mstruct(fwdpar);
else
    error([F ' does not exist in path! STOP.']);
end
F=[name '_Init.mat'];
if exist(F,'file')
    load(F)  ;
    disp([' ']);disp([mfilename,': initial conditions loaded from file ',F]);
    T0=Tinit;
else
    disp(mfilename,': no T initial conditions loaded. equilibrium assumend')
    T0=[];
end

% SITE SPECIFIC SRC PATHS
addpath([locpath]);


% SOLVE FORWARD PROBLEM
dz=diff(z);nz=length(z);
out=0;

if isempty(T0)
    T0=heat1dns(k, kA, kB,r,h,p,qb,POM,dz,ip,maxitnl,tolnl,freeze,out);
end


[Tcalc]=heat1dnt(k,kA,kB,h,r,c,rc,p,qb,...
    dz,ip,dt,it,GST,T0,theta,maxiter,tol,freeze,out)
%==========================================================================
%  Postpocessing
%==========================================================================

% SAVE DATA
filename=strcat([name,'_FwdStat.mat']);
disp(['   ']);disp([' ...save results to: ',filename]);
save(filename,'Tcalc','z')

disp(['   '])

