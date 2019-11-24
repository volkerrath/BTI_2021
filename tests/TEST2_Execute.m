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
plotit = 0;
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

gsth_file           = 'Test2_GSTH.dat';
tlog=2000*y2s;
refyr=tlog;
K      =    2.5;            Ki =  K;       %   = [ 1.  2.  4.];
H      =    0.*1.e-6;      Hi =  H;       %   = [ 0.  2.  4.]*1.e-6;
Qb     =    -60 *1e-3;      Qbi = Qb;      %    = [-40 -60  -80]*1e-3;
nz     =    251;            nzi = nz; % [101 201 301 401 501 1001];
nt     =    251;            nti = nt;      %   = [101 201 301 401 501 1001];


R       = 1000;
C       = 1000;
RC      = R*C;
D       = K/RC;
P       = 0.000001;

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



if run_Nz
    disp([ ' Running test for Nz'])
    TT =  {};
    for nz = nzi

        disp([ newline 'nz = ' num2str(nz)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE Z-MESH
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE MESH OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        F=strcat([name,'_Mesh_in.mat']);
        
        set_z = 1; zstart = 5; zend = 5000; ztype = 'log'; nz=nz;  
        set_t = 1; tstart = 115000*y2s; tend = 10*y2s; ttype= 'log'; nt=nt;
        mesh_in=mstruct(set_z, set_t, ...
            zstart, zend, ztype, nz, tstart, tend, ttype, nt, site);
        save(F,'mesh_in');
        Proc=strcat([site,'_Mesh(name);']); eval(Proc);
        
        F=strcat([name,'_DepthGrid.mat']);
        load(F)
        F=strcat([name,'_TimeGrid.mat']);
        load(F)
                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GSTH SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        gsth_form= 'steps';
        gsth_method = 'linear';
        gsth_file='Test2_GSTH.dat';

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        plotit= 0;
        mod_in=mstruct(plotit,Qb,H,K,gsth0,C,P,R,site);
        F=[name,'_Mod_in'];
        save(F,'mod_in');
        disp(strcat([' generate model for ' name]));
        Proc=strcat([site,'_Mod(name);']); eval(Proc);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE INITIAL VALUES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE INIT OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        plotit=0;
        init_type='equi';
        gsth_form= 'steps';
        gsth_method = 'linear';
        gsth_file='Test2_GSTH.dat';
       
        init_in=mstruct(plotit,init_type,...
            gsth_form,gsth_method,pom,gsth0,gsth_file);
        F=[name,'_Init_in'];
        save(F,'init_in');
        disp(strcat([' generate initial values for ' name]));
        Proc=strcat([site,'_Init(name);']); eval(Proc);
        
              
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature without heat production
        % Beltrami, H. & Mareschal, J.-C.
        %   Recent warming in eastern Canada inferred from Geothermal
        %   Measurements
        %   Geophys. Res. Lett., 1991, 18(4), 605-608
        % Beltrami, H.; Jessop, A. M. & Mareschal, J.-C.
        %   Ground temperature histories in eastern and central Canada from
        %   geothermal measurements: evidence of climate change
        %   Palaeogeography, Palaeoclimatology, Palaeoecology, 1992, 98,
        %   167-184
        % Mareschal, J.-C. & Beltrami, H.
        %   Evidence for recent warming from perturbed thermal gradients:
        %   examples from eastern Canada
        %   Clim. Dyn., 1992, 6, 135-143
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PARAMETER FOR ANALYTICAL FORWARD MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gsth_file='Test2_GSTH.dat';
        gsth_data   =         importdata(gsth_file);
        tgsth       =         gsth_data(:,1)*y2s;
        Tgsth       =         gsth_data(:,2);
        Tgsth       =         [Tgsth; Tgsth(end)];
        pom         =         Tgsth(1)-5.;
        gsth0       =         Tgsth(end);
        gsth_smooth           =         0;
        
        Tgsth(end)
        
        out=0;
        [Tgst]      =         set_stpgst(t,Tgsth,tgsth,gsth_smooth,pom,out);
        
        [Tai,Qai]=heat1dat(K,R,C,Qb,z,t,Tgst,gsth0,tlog,refyr,out);
         Ta = Tai.val;
             
        
        
        
        
        Fwd_tran(name);
%         
        filename=strcat([name,'_FwdTran.mat']);
        disp(['   ']);disp([' ...load results from: ',filename]);
        load(filename)
       

        Tn =Tcalc;
        TT{end+1}=[z(:) Ta(:) Tn(:)];

        F=strcat([name '_N',num2str(nz)]);disp(' ');disp([' Results written to: ', F]);
        save(F, 'Ta','Tn', 'nz');
%         
    end
    
    plotit = 1;
    if plotit
        step=10;
        
        figure    % PLOT TEMPERATURES
        step = 10;
        icurv = 0;
        for nzx=nzi
            icurv=icurv+1;
            Tmp=TT{icurv};
            legstr{icurv}= strcat([' Nz = ',num2str(nzx)]);
            F=strcat([name '_N',num2str(nzx)]);load(F);
            disp([' Results loaded from: ', F]);
            %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
            plot(Tmp(:,2),Tmp(:,1),'-','LineWidth',1); hold on
            plot(Tmp(3*icurv:step:end,3),Tmp(3*icurv:step:end,1),'o','LineWidth',1); hold on
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
        
        figure    % PLOT TEMPERATURE DIFFS
        icurv = 0;
        for nzx =nzi
            icurv=icurv+1;
            Tmp=TT{icurv};
            legstr{icurv}= strcat([' Nz = ',num2str(nzx)]);
            %          strcat([' Nz = ',num2str(nzx),', dx = ',num2str(5000/(nzx-1)),' m']);
            strcat([' Nz = ',num2str(nzx)]);
            F=strcat([name '_N',num2str(nzx)]);load(F);
            Tr=Tmp(:,2)-Tmp(:,3);
            %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
            plot(Tr(:),Tmp(:,1),'-','LineWidth',1); hold on
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