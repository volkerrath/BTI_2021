function [ierr]=OKU_Prep(name)
% Site OTUKUMPU
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


site             = 'ULL';
props           = 'ull';
out             = 0;
estq            = 1;
estq_int        = [100 200 300 400 500];
plotit          = 1;

% INITIAL
theta           =   1.;                   % time steping weight 1/FI .5/CN
maxitnl         =   4;                    % number of nl iterations
tolnl           =   0.00001;
relaxnl         =  1.;
freeze          =  1;                     % include freezing/thawing

% PARAMETER FOR PREPROCESSING
zDatTop         =   75.;
zDatBot         =   1540;

Qb              =  -51e-3;
Qbshift         = -0.0000;
Qb              =   Qb+Qbshift;

GST0            =   7.3;
POM             = -4; 
POR             =   0.01;

Kmin1           =   1.9;
Kmax1           =   4.2;
Kavg1           =   2.6;
rho1            =   2767.;
c1              =   850.;
RC1             =   rho1*c1;

Kmin2           =   3.0;
Kmax2           =   3.7;
Kavg2           =   3.3;
rho2            =   2647.;
c2              =   850.;
RC2             =   rho2*c2;

KTop            =   Kavg1;
KBot            =   Kavg2;

H               =   0.4e-6;
ErrDeflt        =   0.1;

avgmeth ='h';
smooth_data='s';
nspline=65;
wspline=0.6;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE OKU_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_prep_in.mat']);
if exist(F)
    load(F); mstruct(prep_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath([srcpath,filesep,strcat(['src/props/',props])]);


if plotit
    set_graphpars
end

disp(['   ']);
disp(strcat([ ' ...Preprocessing site ', site]));

step = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: SPATIAL AND TEMPORAL MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': set grids']));
% SPATIAL AND TEMPORAL MESH
disp([' ...set grids'])

meshfileZ=[name,'_ZGrid'];
load (meshfileZ);
meshfileT=[name,'_TGrid'];
load (meshfileT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: READ & PREPROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': read obs']));

OT    =   importdata([datpath,'/Ullrigg_March2013.dat']);

zT = OT.data(:,2);zT=(abs(zT));
T = OT.data(:,1);

filename=strcat([datpath site '_obs.mat']);
save(filename, 'T','zT')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),...
    ': average & interpolate to FD mesh']));

% CELL CENTERS
zm=0.5*(z(1:nz-1)+z(2:nz));nc=length(zm);

Tm=NaN(size(zm));
for cell = 1:nc
    incell = find(zT >= z(cell) & zT<z(cell+1));
    if ~isempty(incell)
        Tc=T(incell);
        w=ones(size(Tc));
        Tm(cell) = vavg(Tc,w,'h');
    else
%         disp(strcat(['...>>> cell: ',num2str(cell),' centered at ',...
%                      num2str(zm(cell)),' m empty']))
    end
end
Tm=Tm(isfinite(Tm));

% TEMPERATURE
breaks=linspace(zDatTop,zDatBot,nspline);
Ti = Tm(zm>=zDatTop & zm <=zDatBot);
zi = zm(zm>=zDatTop & zm <=zDatBot);
ss=splinefit(zi,Ti,breaks,'r',wspline);dd=ppdiff(ss);
Tn=ppval(ss,zi);zn=zi;
% 
id=find(z>=zDatTop & z <=zDatBot);
Tobs=Tn;zobs=zn; nd=length(zn);Tobs=Tobs';


KK = Kavg2*ones(size(zm));
KK(zm<=800.) = Kavg1;

RC = RC2*ones(size(zm));
RC(zm<=800.) = RC1;
%

if plotit
    % PLOT TEMPERATURES
    figure
    plot(T,zT,'.b','LineWidth',1); hold on
    plot(Tobs,zobs,'-r','LineWidth',1); hold on
    xlim([0 50]);
    ylim([0 2750]);
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
    
    %
    %
    % PLOT CONDUCTIVITIY
    figure

    plot(KK,zm,'-r','LineWidth',2); hold on
    
    ylim([0 2750]);
    xlim([0 10]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('\lambda (W m^{-1}K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat([site, ' thermal conductivities']),'FontSize',fontsz,'FontWeight',fontwg);
    S1=strcat(['orig']);
    S2=strcat(['mod']);
    legend(S1,S2,'location', 'southeast')
    grid on
    file=strcat([ name '_Lamb']);
    saveas(gcf,file,plotfmt)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP III: ESTIMATE Qb FOR CHOSEN DEPTHINTERVAL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


if estq
    step=step+1;
    disp(strcat([ ' ...>>> Step ',num2str(step),': Estimate Qb - only info']));
    
    zTx=zobs';
    Tx=Tobs;
    TG=diff(Tx)./diff(zTx);
    zG= 0.5*(zTx(1:length(zTx)-1)+zTx(2:length(zTx)));
    Kx=KK(id);Kx=Kx(1:length(Kx)-1)';
    Qobs = Kx.*TG;
    
    
    for Q_EstDepth=zDatBot-estq_int
        DepthInterval=zG>Q_EstDepth;
        Q = Qobs(DepthInterval);
        Q = Q(isfinite(Q));
        Qbm = median(Q);Qbt=mad(Q);
        disp(strcat(['Depth = ', ...
            num2str(Q_EstDepth),' - ',...
            num2str(max(zG)) ' m :', ...
            ' Qmed = ',num2str(Qbm*1000), ...
            ' Qmad = ',num2str(Qbt*1000), ...
            ' mW/m**2' ]));
    end
    
    if plotit
        % PLOT GRADIENTS
        figure
        plot(Qobs*1000,zG,'-r','LineWidth',2); hold on
        ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('HFD T (mW/m^2)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site,' Q_z']),'FontSize',fontsz,'FontWeight',fontwg);
        % S1=strcat(['original obs']);
        %legend(S1,'location', 'southwest')
        grid on
        file=strcat([ name '_Qz']);
        saveas(gcf,file,plotfmt)
        
        figure
        plot(TG*1000,zG,'-r','LineWidth',2); hold on
        ylim([0 2750]);
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        xlabel('\nabla T (mW/m^2)','FontSize',fontsz,'FontWeight',fontwg)
        ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
        title(strcat([site,' Qb']),'FontSize',fontsz,'FontWeight',fontwg);
        % S1=strcat(['original obs']);
        %legend(S1,'location', 'southwest')
        grid on
        file=strcat([ name '_GradT']);
        saveas(gcf,file,plotfmt)
        
    end
end


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
ip=ip;

dz=diff(z);
gts=GST0;
qb=Qb;

k=KK;
% kA=0.0013*nones';
% kB=0.0029*nones'; % OKU fit Kukkonen
kA=0.00*nones';
kB=0.00*nones';

c=0.;r=0.;
rc=RC;  %

h = 0.4*nones';
h(zm>800.)= 1.6e-6;
h(zm>1250.)= 2.7e-6;
h(zm>1340.)= 3.8e-6;
h(zm>zDatBot) = 0.;

p = POR*nones';



% STRUCTURE DATA
nd=length(id);
errT=ErrDeflt;
Terr=errT*ones(size(id))';
Tcov=spdiags(Terr.^2,0,nd,nd);


sitepar=mstruct(k,kA,kB,h,p,c,r,rc,z,ip,t,it,qb,gts,Tobs,id,zobs,Tcov,Terr,props,name);


F=strcat([name '_SITEPar.mat']);
save(F,'sitepar');
disp([' >>>>> site parameter saved to:' F]);
disp([' ']);


end
