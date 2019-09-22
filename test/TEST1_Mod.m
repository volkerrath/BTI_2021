function [ierr]=SITE_Prep(name)
% Prepare Data and Model for inversion

ierr = 0;

load('common.mat')


% SET PATHS

addpath([srcpath,'/src']);
addpath([srcpath,'/tools']);
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
% STEP: SPATIAL MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': get grids']));
meshfileZ=[name,'_DepthGrid.mat'];
load (meshfileZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: PREPROCESS MOEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),...
    ': setup FD mesh']));


% CELL CENTERS
zm=0.5*(z(1:nz-1)+z(2:nz));nc=length(zm);

% BULK THERMAL CONDUCTIVITY, RHOB, RHOC, POR, HEAT PRODUCTION

Kn          = ones(size(dz)).*K;
Pn          = ones(size(dz)).*P; 
Rn          = ones(size(dz)).*R;
Cn          = ones(size(dz)).*C;
Hn          = ones(size(dz)).*H;



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP IV: DEFINE STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': setup structure']));
%


% STRUCTURE MODEL
nz=length(z);
ip=[1:nz-1];nip=length(ip);
nones=ones(nip,1)';
it =1;t=0;

dz=diff(z);
gts=T0;
qb=Qb;

k=K;
kA=0.00*nones';
kB=0.00*nones';


r=R;
c=C;
rc=r.*c;

h = Hn;
p = Pn;


sitemod=mstruct(k,kA,kB,h,p,r,c,rc,z,ip,t,it,qb,gts,props,name);


F=strcat([name '_SiteMod.mat']);
save(F,'sitemod');
disp([' >>>>> site parameter saved to:' F]);
disp([' ']);


end
