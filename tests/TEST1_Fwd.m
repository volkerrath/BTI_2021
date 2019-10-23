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
archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


yeartosec=31557600;sectoyear=1/yeartosec;

run_Q=1;
run_H=1;
run_K=1;
run_Nz=1;



%GRAPHICS
plotit = 1;
set_graphpars
%plotfmt='epsc2';
plotfmt='png';

%
site       = 'TEST1';
props       = 'const';
prepstr       = '_steady_state';

name=[site prepstr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT HALF-SPACE  PHYSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0      =  10.  ; % 10; %

R       = 1000;
C       = 1000;
K       = 2.5;
RC      = R*C;
D       = K/RC;
P       = 0.00000;
Km      = 2.5;          Ki      = [ 1.  2.  4.];
Hm      =  2.*1.e-6;    Hi      = [ 0.  2.  4.]*1.e-6;
Qbm      = -60 *1e-3;   Qbi     = [-40 -60  -80]*1e-3;
nzm     = 251;          nzi     = [101 201 301 401 501 1001];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nz=nzm;
mesh_in=mstruct(set_z, set_t, zstart, zend, ztype, nz, site);
save(F,'mesh_in');
Proc=strcat([site,'_Mesh(name);']);
eval(Proc);

F=strcat([name,'_DepthGrid.mat']);
load(F)

if run_Q 
    disp([ ' Running test for Qb'])
    Ta = [z(:)];Tn = [z(:)];
    for Qb = Qbi
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere, Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zb=zend;
        Tia=T0 - (Qb.*z/Km)  + (Hm.*z/Km).*(zb-z/2);
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
        H = Hm;
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
        disp([' ...load results from: ',filename]);
        load(filename)
        
        Tn =  [Tn Tcalc(:)];
    end
    
    F=strcat([name '_All_Q']);disp(' ');disp([' Results written to: ', F]);
    save(F, 'Ta','Tn', 'Qbi');
    
    if plotit
        step=10;
        
        figure    % PLOT TEMPERATURES
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tn(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' temperatures/Q']),'FontSize',fontsz,'FontWeight',fontwg);
        S1=strcat(['numerical']);
        S2=strcat(['analytic']);
        legend(S1,S2,'location', 'best')
        grid on
        file=strcat([name '_Temp_Q']);
        saveas(gcf,file,plotfmt)
        
        figure    % PLOT TEMPERATURES
        Tr=Tn-Ta;
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tr(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        %plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' residuals/Q']),'FontSize',fontsz,'FontWeight',fontwg);   
        textloc('T_{r} = T_{n}-T_{a}','northeast','Color','blue','FontSize',fontsz,'FontWeight',fontwg);
        grid on
        file=strcat([name '_Res_Q']);
        saveas(gcf,file,plotfmt)
        
    end
end

if run_K
    disp([ ' Running test for K'])
    Ta = [z(:)];Tn = [z(:)];
    for K = Ki
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere, Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zb=zend;
        Tia=T0 - (Qbm.*z/K)  + (Hm.*z/K).*(zb-z/2);
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
        H = Hm;
        Qb = Qbm;
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
        disp([' ...load results from: ',filename]);
        load(filename)
        
        Tn =  [Tn Tcalc(:)];
    end
    
    F=strcat([name '_All_K']);disp(' ');disp([' Results written to: ', F]);
    save(F, 'Ta','Tn', 'Ki');
    
    if plotit
        step=10;
        
        figure    % PLOT TEMPERATURES
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tn(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' temperatures/Q']),'FontSize',fontsz,'FontWeight',fontwg);
        S1=strcat(['numerical']);
        S2=strcat(['analytic']);
        legend(S1,S2,'location', 'best')
        grid on
        file=strcat([name '_Temp_K']);
        saveas(gcf,file,plotfmt)
        
        figure    % PLOT TEMPERATURES
        Tr=Tn-Ta;
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tr(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        %plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' residuals/Q']),'FontSize',fontsz,'FontWeight',fontwg);
        textloc('T_{r} = T_{n}-T_{a}','northeast','Color','blue','FontSize',fontsz,'FontWeight',fontwg);
        grid on
        file=strcat([name '_Res_K']);
        saveas(gcf,file,plotfmt)
        
    end
end

if run_H
    disp([ ' Running test for H'])
    Ta = [z(:)];Tn = [z(:)];
    for H = Hi
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere, Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zb=zend;
        Tia=T0 - (Qbm.*z/Km)  + (H.*z/Km).*(zb-z/2);
        disp([ 'H_{in} = ' num2str(H*1e6)])
        Ta = [Ta Tia(:)];
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL FORWARD MODEL SETUP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % VARIABLES SET HERE OUTSIDE PREP OVERWRITE DEFAULTS INSIDE!
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Qb = Qbm;
        K  = Km;
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
    end
    
    F=strcat([name '_All_H']);disp(' ');disp([' Results written to: ', F]);
    save(F, 'Ta','Tn', 'Hi');
    
    if plotit
        step=10;
        
        figure    % PLOT TEMPERATURES
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tn(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' temperatures/H']),'FontSize',fontsz,'FontWeight',fontwg);
        S1=strcat(['numerical']);
        S2=strcat(['analytic']);
        legend(S1,S2,'location', 'best')
        grid on
        file=strcat([name '_Temp_H']);
        saveas(gcf,file,plotfmt)
        
        figure    % PLOT TEMPERATURES
        Tr=Tn-Ta;
        %plot(Ta(:,2:end),Ta(:,1),':','LineWidth',1); hold on
        plot(Tr(:,2:end),Tn(:,1),'-','LineWidth',1); hold on
        %plot(Ta(1:step:end,2:end),Ta(1:step:end,1),'o','LineWidth',1); hold on
        %     xlim([0 50]);
        %     ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site, ' residuals/H']),'FontSize',fontsz,'FontWeight',fontwg);
        textloc('T_{r} = T_{n}-T_{a}','northeast','Color','blue','FontSize',fontsz,'FontWeight',fontwg);
        grid on
        file=strcat([name '_Res_H']);
        saveas(gcf,file,plotfmt)
        
    end
end

if run_Nz
    disp([ ' Running test for Nz'])
    for nz = nzi
        
        disp([ 'nz_{in} = ' num2str(nz)])
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
        mesh_in=mstruct(set_z, set_t, zstart, zend, ztype, nz, site);
        save(F,'mesh_in');
        Proc=strcat([site,'_Mesh(name);']);
        eval(Proc);
        
        F=strcat([name,'_DepthGrid.mat']);
        load(F)
        Ta = [z(:)];Tn = [z(:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL FORWARD MODEL
        % Analytical  solution for temperature with heat production
        % Stüwe, K. Geodynamiics of the Lithosphere, Springer, Berlin, 2002
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zb=zend;
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