
clear all
close all
clc

% SET PATHS
pltpath='./';
datpath='./';
srcpath='./';
utlpath='./';

save('common','srcpath','utlpath','datpath','pltpath'),

if ~exist(pltpath,'dir'),
    mkdir(pltpath)
end

addpath([pwd,filesep,'./src']);
addpath([pwd,filesep,strcat(['./tools'])]);


dfmt=1;ffmt='.zip';
archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);
yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS
%plotfmt='epsc2';
plotfmt='png';

linwdt=2;
fontwg='normal';

fontsz=16;
plot_log=1;
plot_gst=0;
plot_prep=0;
plot_btp=0;
plot_btq=0;
plot_btr=0;
plot_btg=0;
pltdepth       =   1500;


run_parallel=1;
parcors=   4;

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

%
sites       = {'LX2','LAX'};
props       = {'lax','olk'};
%nsites=length(sites);
whichsites=[1 2];%


filename=strcat([sites{whichsites},'_common.mat'])
save(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FOR FWD CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['   ']);
disp([' ...set up meshes for paleoclimate inversion  ']);


meshfile='LAX_mesh';
load (meshfile)
nz = length(z);
nt = length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE COMMON NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta           =   1.*ones(1,nt);              % time steping weight 1/FI .5/CN
%     theta          = 0.5*ones(1,nt);
%     theta(1:10)    =  1.;
% NONLINEAR ITERATION
maxitnl         =   4;                          % number of nl iterations
tolnl           =   0.00001;
relaxnl         =  1.;
freeze          =  1;                           % include freezing/thawing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FOR INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' ']);disp([' ...set up parametrization for paleoclimate inversion  ']);

% DEFINE LOGARITHMIC INVERSION GRID
%IPOM=-4;
nsteps          =   25;                 % number of steps
base            =   0.;                 % base
tstart          =  115000*yeartosec;     % from
tend            =  10*yeartosec;         % to
[pt,it]=set_mgsth(t,base,tstart,tend,nsteps);

% SENSITVITY CALCULATION PARAMETERS
diffmeth        =   'FD';               % how to calculate Sensitivity
dp              =   0.001;               % fd PERTURBATION

% ITERATION PARAMETER
tol_solve       =   0.00001;             % tolerance for solver
maxiter_solve   =   32;                 % number of solver iterations
tol_inv         =   [0.0001,0.00001];  % tolerance for inversion rms
maxiter_inv     =   100;     % maximal number of iterations




% % REGULARIZATION PARAMETER FOR TIKHONOV
% reg_opt          =  'UPRE';             % Unbiased Predictive Risk Estimator
% reg_opt           =   'GCV';          % Generalized Cross-Validation
% reg_opt           =   'LC';               % L-Curve
reg_opt         =   'FIX';              % fixed regularization parameter
start_regpar    =    maxiter_inv+1;                 % when to start search
modul_regpar    =    1;                 % how often
mregpar_adaptint =    1;
msteps_regpar    =   48;                % number of test values

regpar0=[1  1 1];
reg0par=[0.1 0.01  0]; %logspace(-1,1,10);31
reg1par=[0.01];  %logspace(-3.,2,msteps_regpar);
reg2par=[0];
% RELAXATION
relax           =   1;% 0.98;
start_relax     =   2;
modul_relax     =   1;
min_relax       =   0.1;
reg_shift=1;
outsteps=0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE/LOAD  SITE-SPECIFIC SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


smooth_log=0;

Q0=-[46]*1e-3;


smooth_props='m';
avgmeth ='h';
smooth_log=0;
smooth_data='f';
LFD=7;TFD={'tri','mir'};MD=1;
smooth_grads = 'm';

prepstr = '';
start_regpar    =   0;


regopts={'GCV'};
reg1par=logspace(-2,3,msteps_regpar);
for iprop = 1:length(props)
    prop_opt=props{iprop};PROP=prop_opt;
    addpath([pwd,filesep,strcat(['./src/props/',PROP])]);
    for ireg=1:length(regopts)
        reg_opt=regopts{ireg};
        for reg0par=[ 0.003 ];
            for H = [2.e-6 1.4e-6 ]
                for reg2par=[0];
                    for nsteps=[24];
                        [pt,it]=set_mgsth(t,base,tstart,tend,nsteps);
                        for m_apr_set=[4];
                            m_ini_set=[m_apr_set];
                            for isites=whichsites
                                SITE=sites{isites};
                                NAME=strcat([SITE,...
                                    '_Tikh',reg_opt,...
                                    '_Prior',num2str(m_apr_set),...
                                    '_H',num2str(H*1e7),...
                                    '_Props_',PROP]);
                                
                                
                                Tinitial= strcat(['LAX-',PROP,'_H',num2str(H*1e7),'_initial']);
                                % INITIAL
                                if exist('Tinitial')
                                    load(Tinitial);
                                    disp([' ']);disp([' initial condition loaded from file ',Tinitial]);
                                    T0=Tinit;
                                else
                                    T0=[];
                                end
                                numpar=mstruct(T0,theta,maxitnl,tolnl,relaxnl,freeze);
                                %regpar=setregpar(regpar0,reg0par,reg1par,reg2par);
                                
                                disp([' ']);disp([' ...set up site-specific parameter settings  ']);
                                filename=strcat([NAME,'_invpar.mat']);
                                save (filename);
                                
                                S=sites{whichsites(isites)}; N=NAME;P=props{isites};
                                disp([' ']);disp(strcat([' generate model for ' N ]));
                                C=strcat([ S ,'_Prep',prepstr,'(S,N,P);']); eval(C);
                                disp([' ']);
                                
                                
                                
                                % GSTH_HybrS(SITE,NAME);
                                GSTH_TikhSX(SITE,NAME,PROP);
                                % rmpath([srcpath,strcat(['src/props/',PROP])]);
                                %close all
                            end
                        end
                    end
                end
            end
        end
    end
    rmpath([pwd,filesep,strcat(['./src/props/',PROP])]);
end


regopts={'FIX'};
H = 2.e-6;
for iprop = 1:length(props)
    prop_opt=props{iprop};PROP=prop_opt;
    addpath([pwd,filesep,strcat(['./src/props/',PROP])]);
    for ireg=1:length(regopts)
        reg_opt=regopts{ireg};
        for reg0par=[ 0.003 ];
            for reg1par=[3 0.3 0.03 0.003];
                for reg2par=[0];
                    for nsteps=[24];
                        [pt,it]=set_mgsth(t,base,tstart,tend,nsteps);
                        for m_apr_set=[6 4 2];
                            m_ini_set=[m_apr_set];
                            for isites=whichsites
                                SITE=sites{isites};
                                NAME=strcat([SITE,...
                                    '_Tikh',reg_opt,...
                                    '_RegP_',num2str(reg0par,'%4.3f'),...
                                    '_',num2str(reg1par,'%4.3f'),...
                                    '_Prior',num2str(m_apr_set),...
                                    '_H',num2str(H*1e7),...
                                    '_Props_',PROP]);
                                
                                
                                Tinitial= strcat(['LAX-',PROP,'_H',num2str(H*1e7),'_initial']);
                                % INITIAL
                                if exist('Tinitial')
                                    load(Tinitial);
                                    disp([' ']);disp([' initial condition loaded from file ',Tinitial]);
                                    T0=Tinit;
                                else
                                    T0=[];
                                end
                                numpar=mstruct(T0,theta,maxitnl,tolnl,relaxnl,freeze);
                                %regpar=setregpar(regpar0,reg0par,reg1par,reg2par);
                                
                                disp([' ']);disp([' ...set up site-specific parameter settings  ']);
                                filename=strcat([NAME,'_invpar.mat']);
                                save (filename);
                                
                                S=sites{whichsites(isites)}; N=NAME;P=props{isites};
                                disp([' ']);disp(strcat([' generate model for ' N ]));
                                C=strcat([ S ,'_Prep',prepstr,'(S,N,P);']); eval(C);
                                disp([' ']);
                                
                                
                                
                                % GSTH_HybrS(SITE,NAME);
                                GSTH_TikhSX(SITE,NAME,PROP);
                                % rmpath([srcpath,strcat(['src/props/',PROP])]);
                                %close all
                            end
                        end
                    end
                end
            end
        end
    end
    rmpath([pwd,filesep,strcat(['./src/props/',PROP])]); 
end
% H = 1.4e-6;

% for iprop = 1:length(props)
%     prop_opt=props{iprop};PROP=prop_opt;
%     addpath([pwd,filesep,strcat(['./src/props/',PROP])]);
%     for ireg=1:length(regopts)
%         reg_opt=regopts{ireg};
%         for reg0par=[ 0.003 ];
%             for reg1par=[3 0.3 0.03 0.003];
%                 for reg2par=[0];
%                     for nsteps=[24];
%                         [pt,it]=set_mgsth(t,base,tstart,tend,nsteps);
%                         for m_apr_set=[6 4 2];
%                             m_ini_set=[m_apr_set];
%                             for isites=whichsites
%                                 SITE=sites{isites};
%                                 NAME=strcat([SITE,...
%                                     '_Tikh',reg_opt,...
%                                     '_RegP_',num2str(reg0par,'%4.3f'),...
%                                     '_',num2str(reg1par,'%4.3f'),...
%                                     '_Prior',num2str(m_apr_set),...
%                                     '_H',num2str(H*1e7),...
%                                     '_Props_',PROP]);
%                                 
%                                 
%                                 Tinitial= strcat(['LAX-',PROP,'_H',num2str(H*1e7),'_initial']);
%                                 % INITIAL
%                                 if exist('Tinitial')
%                                     load(Tinitial);
%                                     disp([' ']);disp([' initial condition loaded from file ',Tinitial]);
%                                     T0=Tinit;
%                                 else
%                                     T0=[];
%                                 end
%                                 numpar=mstruct(T0,theta,maxitnl,tolnl,relaxnl,freeze);
%                                 %regpar=setregpar(regpar0,reg0par,reg1par,reg2par);
%                                 
%                                 disp([' ']);disp([' ...set up site-specific parameter settings  ']);
%                                 filename=strcat([NAME,'_invpar.mat']);
%                                 save (filename);
%                                 
%                                 S=sites{whichsites(isites)}; N=NAME;P=props{isites};
%                                 disp([' ']);disp(strcat([' generate model for ' N ]));
%                                 C=strcat([ S ,'_Prep',prepstr,'(S,N,P);']); eval(C);
%                                 disp([' ']);
%                                 
%                                 
%                                 
%                                 % GSTH_HybrS(SITE,NAME);
%                                 GSTH_TikhSX(SITE,NAME,PROP);
%                                 % rmpath([srcpath,strcat(['src/props/',PROP])]);
%                                 %close all
%                                 
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%    rmpath([pwd,filesep,strcat(['./src/props/',PROP])]);
% end
