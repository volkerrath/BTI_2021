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

%GRAPHICS

set_graphpars
%plotfmt='epsc2';
plotfmt='png';
plotit = 1;
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
Hm      =  2.*1.e-6;
Hi      = [ 0.  2.  4.]*1.e-6;
Qbm      = -60 *1e-3;
Qbi     = [-40 -60  -80]*1e-3;
nzm     = 251;
nzi     = [101 201 301 401 501 1001];
ntm     = 301;
nti     = [201 301 401 501 601 1001];


for nz = nzi
    
    disp([ 'nz_{in} = ' num2str(nz)])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERATE JOINT Z-MESH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % VARIABLES SET HERE OUTSIDE MESH OVERWRITE DEFAULTS INSIDE!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    F=strcat([name,'_Mesh_in.mat']);
    set_z   = 1;
    zstart  = 0;
    zend    = 5000;
    ztype   = 'linear';
    set_t   = 1;
    tstart  = 3;
    tzend    = 120000;
    ttype   = 'log';
    
    mesh_in=mstruct(...
        set_z, zstart, zend, ztype, nz,...
        set_t, tstart, tend, ttype, nt, site);
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
    Tia=T0 - (Qbm.*z/K)  + (Hm.*z/K).*(zb-z/2);
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
    file=strcat([name '_TempN']);
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
    legend(legstr,'location','best');
    grid on
    file=strcat([name '_Res_N']);
    saveas(gcf,file,plotfmt)
    
end

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
    Tia=T0 - (Qbm.*z/K)  + (Hm.*z/K).*(zb-z/2);
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
    for ntx=nti         
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
    title(strcat([site, ' temperatures/Nt']),'FontSize',fontsz,'FontWeight',fontwg);
    legend(legstr,'location','best');
    grid on
    file=strcat([name '_TempNt']);
    saveas(gcf,file,plotfmt)
    
    figure    % PLOT TEMPERATURES
    icurv = 0;
    for ntx =nzi
         icurv=icurv+1;
         legstr{icurv}= strcat([' Nt = ',num2str(ntx)]);
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
    legend(legstr,'location','best');
    grid on
    file=strcat([name '_Res_Nt']);
    saveas(gcf,file,plotfmt)
    
end
