clear all
close all
clc

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

% SET PATHS
pltpath='./';
datpath='./';
srcpath='../../../';
utlpath='../../../';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

% ONLY FOR PARRALLEL  EXECUTION
run_parallel=1;
parcors=   8;



save('common','srcpath','utlpath','datpath','pltpath','run_parallel','parcors'),



dfmt=1;ffmt='.zip';
%archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);
yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS
plot_results = 0;

set_graphpars
%plotfmt='epsc2';
plotfmt='png';

%
site       = 'SYNA';
props       = 'syn';


Qb = -30*2.3253*1e-3;

NSamp=32;
for sample = [1:NSamp]
    for ErrDeflt = [0.03 0.1 0.3 1.]
        
        prepstr       = strcat( ['_Sample',num2str(sample),'_Err',num2str(ErrDeflt*1000,'%.2g'),'mK'] );
        
        name=[site prepstr];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL PARAMETER FOR FORWARD MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta           =   1.;                   % time steping weight 1/FI .5/CN
        maxitnl         =   4;                    % number of nl iterations
        tolnl           =   0.00001;
        relaxnl         =  1.;
        freeze          =  1;                     % include freezing/thawing
        
        fwdpar=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
        F=strcat([name,'_FwdPar.mat']);
        save(F, 'fwdpar')
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE MESHES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE ULL_MESH OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        F=strcat([name,'_Mesh_in.mat']);
        set_z = 1;
        set_t = 1;
        mesh_in=mstruct(set_z, set_t);
        save(F,'mesh_in');
        disp(strcat([' generate meshes for ' name]));
        C=strcat([site,'_Mesh(name);']);
        eval(C);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE PHYSICAL MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        plotit=0;
        CovType='g';
        %ErrDeflt=0.05;
        L=3;
        zDatTop=10;zDatBot=1500;
        prep_in=mstruct(plotit,CovType,ErrDeflt,L,zDatTop,zDatBot,Qb);
        
        F=[name,'_Prep_in'];
        save(F,'prep_in');
        disp(strcat([' generate model for ' name]));
        C=strcat([site,'_Prep(name);']);
        eval(C);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE INITIAL VALUES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE INIT OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        plotit=0;
        init_type='p';
        init_form= 'steps';
        method = 'linear';
        GSTH_file='GSTHBallingA.csv';
        init_in=mstruct(plotit,init_type,init_form,method,GSTH_file);
        F=[name,'_Init_in'];
        save(F,'init_in');
        disp(strcat([' generate initial values for ' name]));
        C=strcat([site,'_Init(name);']);eval(C);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE INVERSION PARAMETER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F=[name,'_InvPar_in'];
        invpar_in=struct();
        save(F,'invpar_in');
        disp(strcat([' generate inversion setup for ' name]));
        C=strcat([site,'_InvPar(name);']);eval(C);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RUN inversion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        Tikh_gsth(name);
        
        
    end
end
if plot_results
    SYNB_TikhPlot
end
