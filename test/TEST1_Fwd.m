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
plotit = 1;
%
site       = 'TEST1';
props       = 'syn';
prepstr       = '_steadt_state';
name=[site prepstr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT HAKF-SPACE  PHYSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0      =  10.  ; % 10; %
QBi      = [-40 -60  -80]*1e-3;
R       = 1000;
C       = 1000;
K       = 2.5;
RC      = R*C;
D       = K/RC;
P       = 0.00000;
Hi       =  [ 0.  2.  4.]*1.e-6;

name = strcat([ site '_steady-state']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL FORWARD MODEL CTRL PARAMETER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta           =   1.;                   % time steping weight 1/FI .5/CN
maxitnl         =   4;                    % number of nl iterations
tolnl           =   0.00001;
relaxnl         =  1.;
freeze          =  1;                     % include freezing/thawing

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
nz=251;
mesh_in=mstruct(set_z, set_t, zstart, zend, ztype, nz, site);
save(F,'mesh_in');
Proc=strcat([site,'_Mesh(name);']);
eval(Proc);

F=strcat([name,'_DepthGrid.mat']);
load(F)
Ta = [z(:)];Tn = [z(:)];

for Qb = QBi
    for H = Hi
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zb=zend;
        Tia=T0 - (Qb.*z/K)  + (H.*z/K).*(zb-z/2);
        Qia = K*diff(Tia)./diff(z);
        disp([ 'Qb_{in} = ' num2str(Qia(nz-1)*1000)])
        Ta = [Ta Tia(:)];
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mod_in=mstruct(Qb,H,K,T0,C,P,R,site);
        F=[name,'_Mod_in'];
        save(F,'mod_in');
        disp(strcat([' generate model for ' name]));
        Proc=strcat([site,'_Mod(name);']);
        eval(Proc);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Fwd_stat(name);
        Tn =  Ta; %[Tn Tin(:)]
    end
end

F=strcat([name '_All']);disp(' ');disp([' Results written to: ', F]);
save(F, 'Ta','Tn', 'QB','H');

if plotit
    step=50;
    figure    % PLOT TEMPERATURES
    plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
    plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
    %     xlim([0 50]);
    %     ylim([0 2750]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat([site, ' temperatures']),'FontSize',fontsz,'FontWeight',fontwg);
    S1=strcat(['orig']);
    S2=strcat(['mod']);
    legend(S1,S2,'location', 'southwest')
    grid on
    file=strcat([name '_Temp']);
    saveas(gcf,file,plotfmt)
    
    
end