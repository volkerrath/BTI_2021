clear all
close all
clc

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

% SET PATHS
pltpath='./';
datpath='./';
srcpath='../';
utlpath='../';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

% ONLY FOR PARRALLEL  EXECUTION
run_parallel=1;
parcors=   2;


save('common','srcpath','utlpath','datpath','pltpath','run_parallel','parcors'),


dfmt=1;ffmt='.zip';
%archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS

set_graphpars
%plotfmt='epsc2';
plotfmt='png';

%
site       = 'TEST1';
props       = 'syn';
prepstr       = '_steadt_state';
name=[site prepstr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT HAKF-SPACE  PHYSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0      =  10.  ; % 10; %
QB      = [-40 -60  -80]*1e-3;
R       = 1000;
C       = 1000;
K       = 2.5;
RC      = R*C;
D       = K/RC;
H       =  [ 0.  2.  4.]*1.e-6;

name = strcat([ site '_steady-state']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE JOINT Z-MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES SET HERE OUTSIDE MESH OVERWRITE DEFAULTS INSIDE!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F=strcat([name,'_Mesh_in.mat']);
set_z   = 1;set_t   = 0;
zstart  = 0;
zend    = 5000;
ztype   = 'linear';
mesh_in=mstruct(set_z, set_t, zstart, zend, ztype);
save(F,'mesh_in');
disp(strcat([' generate meshes for ' name]));
C=strcat([site,'_Mesh(name);']);
eval(C);

F=strcat([name,'_DepthGrid.mat']);
load(F)
Ta = [z(:)];Tn = [z(:)];

for Qi = QB
    for Hi = H
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=T0+(Qi*z/k)  + (H*z/lambda).*(zb-z/2);

        
        Ta = [Ta Tia(:)]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta           =   1.;                   % time steping weight 1/FI .5/CN
        maxitnl         =   4;                    % number of nl iterations
        tolnl           =   0.00001;
        relaxnl         =  1.;
        freeze          =  1;                     % include freezing/thawing
        
        fwdpar=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
        F=strcat([name,'_FwdPar.mat']);
        save(F, 'fwdpar')
        %
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         % GENERATE PHYSICAL MODEL
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %         % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %         %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %         plotit=0;
        %         CovType='g';
        %         ErrDeflt=0.00;
        %         L=0;
        %         zDatTop=0;zDatBot=5000;
        %         prep_in=mstruct(plotit,CovType,ErrDeflt,L,zDatTop,zDatBot,Qb);
        %
        %         F=[name,'_Prep_in'];
        %         save(F,'prep_in');
        %         disp(strcat([' generate model for ' name]));
        %         C=strcat([site,'_Prep(name);']);
        %         eval(C);
        %
        %
        %         Fwd_gsth(name);
                 Tn =  Ta; %[Tn Tin(:)]
    end
end

F=strcat([name '_All']);
save(F, 'Ta','Tn', 'QB','H');