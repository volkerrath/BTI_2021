clear all
close all
clc

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

y2s=3600*24*365.25;s2y=1./y2s;

% SET PATHS
pltpath='./';
datpath='./';
srcpath='../';
utlpath='../';
locpath='./local';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

% ONLY FOR PARRALLEL  EXECUTION
run_parallel=1;
parcors=   2;


save('common','srcpath','utlpath','datpath','pltpath','locpath',...
    'run_parallel','parcors'),


dfmt=1;ffmt='.zip';
% archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


yeartosec=31557600;sectoyear=1/yeartosec;


run_Nt=1;
run_Nz=1;
run_Sp=1;

%GRAPHICS
plotit = 1;
set_graphpars
%plotfmt='epsc2';
plotfmt='png';

%
site       = 'TEST2';
prepstr       = '_transient';

name=[site prepstr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT HALF-SPACE  PHYSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0      =  10.  ; % 10; %

Km      = 2.5;          Ki =  Km;       %   = [ 1.  2.  4.];
Hm      =  2.*1.e-6;    Hi =  Hm;       %   = [ 0.  2.  4.]*1.e-6;
Qbm      = -60 *1e-3;   Qbi = Qbm;      %    = [-40 -60  -80]*1e-3;
nzm     = 251;          nzi = nzm;      %   = [101 201 301 401 501 1001];
ntm     = 251;          nti = ntm;      %   = [101 201 301 401 501 1001];
GSTH_file           = 'Test2_GSTH.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL FORWARD MODEL CTRL PARAMETER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta           =   1.;                   % time steping weight 1/FI .5/CN
maxitnl         =   4;                    % number of nl iterations
tolnl           =   0.00001;
relaxnl         =  1.;
freeze          =  1;                     % include freezing/thawing
fwdpar=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
F=strcat([name,'_FwdPar.mat']);
save(F, 'fwdpar')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE JOINT Z-MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES SET HERE OUTSIDE MESH OVERWRITE DEFAULTS INSIDE!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F=strcat([name,'_Mesh_in.mat']);
set_z   = 1;set_t   = 1;


zstart  = 0;
zend    = 5000;
ztype   = 'log';
nz=nzm;

tstart  = 110000*y2s;
tend    = 30*y2s;
ttype= 'log';
nt=ntm;

mesh_in=mstruct(set_z, set_t, ...
    zstart, zend, ztype, nz,...
    tstart, tend, ttype, nt, site);
save(F,'mesh_in');
Proc=strcat([site,'_Mesh(name);']);
eval(Proc);

F=strcat([name,'_DepthGrid.mat']);
load(F)
F=strcat([name,'_TimeGrid.mat']);
load(F)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FOR ANALYTICAL FORWARD MODEL
% Analytical  solution for temperature without heat production
% Beltrami, H. & Mareschal, J.-C. 
%   Recent warming in eastern Canada inferred from Geothermal Measurements 
%   Geophys. Res. Lett., 1991, 18(4), 605-608
% Beltrami, H.; Jessop, A. M. & Mareschal, J.-C. 
%   Ground temperature histories in eastern and central Canada from 
%   geothermal measurements: evidence of climate change 
%   Palaeogeography, Palaeoclimatology, Palaeoecology, 1992, 98, 167-184
% Mareschal, J.-C. & Beltrami, H. 
%   Evidence for recent warming from perturbed thermal gradients: examples 
%   from eastern Canada Clim. Dyn., 1992, 6, 135-143
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GSTH_data   =         importdata(GSTH_file);
tGSTH       =         GSTH_data(:,1)*y2s;
TGSTH       =         GSTH_data(:,2);
TGSTH       =         [TGSTH; TGSTH(end)]; 
POM         =         TGSTH(1)-4;
T0          =         TGSTH(end);
[Tgst]      =         set_stpgst(t,TGSTH,tGSTH,L,POM,0);


 
if run_Nz
    disp([ ' Running test for Nz'])
    Ta =  [];
    for nz = nzi
        
        disp([ 'nz_{in} = ' num2str(nz)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE Z-MESH
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE MESH OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        F=strcat([name,'_Mesh_in.mat']);
        
        set_z   = 1;set_t   = 1;
                
        zstart  = 0;
        zend    = 5000;
        ztype   = 'log';
        nz=nzi;
        
        tstart  = 110000*y2s;
        tend    = 30*y2s;
        ttype= 'log';
        nt=ntm;
        
        mesh_in=mstruct(set_z, set_t, ...
            zstart, zend, ztype, nz,...
            tstart, tend, ttype, nt, site);
        
        save(F,'mesh_in');
        Proc=strcat([site,'_Mesh(name);']);
        eval(Proc);
        
        F=strcat([name,'_DepthGrid.mat']);
        load(F)
        F=strcat([name,'_TimeGrid.mat']);
        load(F)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere, Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T,Q]=heat1dat(Km,R,C,Qbm,z,t,Tgst,T0,tlog,refyr,out)b=zend;
        Tia=T0 - (Qbm.*z/Km)  + (Hm.*z/Km).*(zb-z/2);
        Ta =  [Ta Tia(:)];
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Qb = Qbm;
        H  = Hm;
        K = Km;
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
        
        filename=strcat([name,'_FwdStat.mat']);
        disp(['   ']);disp([' ...load results from: ',filename]);
        load(filename)
        
        Tn =  [Tn Tcalc(:)];
        F=strcat([name '_N',num2str(nz)]);disp(' ');disp([' Results written to: ', F]);
        save(F, 'Ta','Tn', 'nz');
        
    end
    
    
    
    if plotit
        step=10;
        
        figure    % PLOT TEMPERATURES
        step = 10;
        icurv = 0;
        for nzx=nzi
            icurv=icurv+1;
            legstr{icurv}= strcat([' Nz = ',num2str(nzx)]);
            F=strcat([name '_N',num2str(nzx)]);load(F);
            disp([' Results loaded from: ', F]);
            %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
            plot(Tn(:,2),Tn(:,1),'-','LineWidth',1); hold on
            plot(Ta(1:step:end,2),Ta(1:step:end,1),'o','LineWidth',1); hold on
        end
        %     xlim([0 50]);
        %     ylim([0 2750]);
        
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' temperatures/Nz']),'FontSize',fontsz,'FontWeight',fontwg);
        legend(legstr,'location','best');
        grid on
        file=strcat([name '_Temp_Nz']);
        saveas(gcf,file,plotfmt)
        
        figure    % PLOT TEMPERATURES
        icurv = 0;
        for nzx =nzi
            icurv=icurv+1;
            legstr{icurv}= strcat([' Nz = ',num2str(nzx)]);
            %          strcat([' Nz = ',num2str(nzx),', dx = ',num2str(5000/(nzx-1)),' m']);
            strcat([' Nz = ',num2str(nzx)]);
            F=strcat([name '_N',num2str(nzx)]);load(F);
            Tr=Tn(:,2)-Ta(:,2);
            %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
            plot(Tr(:),Tn(:,1),'-','LineWidth',1); hold on
            %plot(Ta(1:step:end,2),Ta(1:step:end,1),'o','LineWidth',1); hold on
        end
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' residuals/Nz']),'FontSize',fontsz,'FontWeight',fontwg);
        textloc('T_{r} = T_{n}-T_{a}','northeast','Color','blue','FontSize',fontsz,'FontWeight',fontwg);
        legend(legstr,'location','best');
        grid on
        file=strcat([name '_Res_Nz']);
        saveas(gcf,file,plotfmt)
        
    end
end