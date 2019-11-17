function [ierr]=SITE_Mod(name)
% Prepare Data and Model for inversion


ierr = 0;

load('common.mat')

y2s=3600*24*365.25;s2y=1./y2s;
% SET PATHS

addpath([srcpath,'/src']);
addpath([srcpath,'/tools']);
addpath([strcat(['./local'])]);
addpath([datpath]);
% 
% dfmt=1;ffmt='.zip';
% archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);

% GENERAL SETTINGS
out             = 0;
plotit          = 0;

Qb              =  -60*1e-3;
H               =   0.;
T0              =   10.;
P               =   0.00001;
C               =  1000.;
R               =  1000.;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE SYN_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_Mod_in.mat']);
if exist(F)
    disp([' ...',mfilename ' defaults overwritten from ', F])
    load(F); mstruct(mod_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if plotit
    set_graphpars
end

disp(['   ']);
disp(strcat([ ' ...Preprocessing site ', site]));

step = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: MESHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': get grids']));
F=[name,'_DepthGrid.mat'];
load (F);
F=[name,'_TimeGrid.mat'];
load (F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: PREPROCESS MOEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),...
    ': setup FD mesh']));


% CELL CENTERS
zm=0.5*(z(1:nz-1)+z(2:nz));nc=length(zm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: DEFINE SITEMOD STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gsth_form           = 'steps';
gsth_method         = 'linear';
gsth_file           = 'Test2_GSTH.dat';
gsth_smooth         = 0;

gsth_data    =        importdata(gsth_file);
tgsth       =         gsth_data(:,1)*y2s;
Tgsth       =         gsth_data(:,2);
Tgsth       =         [Tgsth; Tgsth(end)];
pom         =         Tgsth(1)-5.;
[Tgst] = set_stpgst(t,Tgsth,tgsth,gsth_smooth,pom,0);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: DEFINE SITEMOD STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': setup sitemod structure']));
%


% STRUCTURE MODEL


qb=Qb;

ip=[1:nz-1];
it=[1:nt-1];
nones=ones(nz,1);nones=nones(:);

k           = nones.*K;
kA          = 0.00*nones;
kB          = 0.00*nones;
r           = nones.*R;
c           = nones.*C;
h           = nones.*H;
p           = nones.*P;
rc          = r.*c;

sitemod=mstruct(k,kA,kB,h,p,c,r,rc,z,ip,t,it,qb,Tgst,pom,name);


F=strcat([name '_SiteMod.mat']);
save(F,'sitemod');
disp([' >>>>> site parameter saved to:' F]);
disp([' ']);


end
