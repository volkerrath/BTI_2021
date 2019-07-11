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


site             = 'OKU';
props           = 'oku';
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
zDatTop         =   100.;
zDatBot         =   1800;

Qb              =  -42e-3;
Qbshift         = -0.0000;
Qb              =   Qb+Qbshift;

GST0            =   0;
POM             = -4; 
POR             =   0.01;
KTop            =   2.3;
KBot            =   2.9;
KMean           =   2.7;
ErrDeflt        =   0.1;

smooth_props='m';
avgmeth ='h';
smooth_data='s';
nspline=101;
wspline=0.5;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE OKU_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_Prep_in.mat']);
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

meshfileZ=[name,'_DepthGrid'];
load (meshfileZ);
meshfileT=[name,'_TimeGrid'];
load (meshfileT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: READ & PREPROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': read obs']));

OT    =   importdata([datpath,'/OKU_Temp_orig.dat']);
OK    =   importdata([datpath,'/OKU_Lamb_orig.dat']);

if isstruct(OT)
    % Version 2012b
    zT = OT.obs(:,1);
    T = OT.obs(:,2);
    zK = OK.obs(:,1);
    K = OK.obs(:,2);
else
    % Version 2010b
    zT = OT(:,1);
    T = OT(:,2);
    zK = OK(:,1);
    K = OK(:,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),...
    ': average & interpolate to FD mesh']));


zKBot=max(max(zK));
zKTop=min(min(zK));

% CELL CENTERS
zm=0.5*(z(1:nz-1)+z(2:nz));nc=length(zm);

% TEMPERATURE
breaks=linspace(min(zT),max(zT),nspline);
ss=splinefit(zT,T,breaks,'r',wspline);dd=ppdiff(ss);
Tx=ppval(ss,z);Tc=ppval(ss,zm);dTc=ppval(dd,zm);
resT=ppval(ss,zT)-T;no=length(resT);errT=norm(resT)/sqrt(no);
Ts=NaN(size(z));
Ts(z>=min(zT) & z <=max(zT))=Tx(z>=min(zT) & z <=max(zT));


% CUT LOG
Ts(z< zDatTop)=NaN;
Ts(z> zDatBot)=NaN;
id=find(isfinite(Ts));
Tobs=Ts(id);zobs=z(id); nd=length(id);Tobs=Tobs';


% BULK THERMAL CONDUCTIVITY
iKok=K>=0.25 & K<=15; K = K(iKok); zK = zK(iKok);
Ks=NaN(size(zm));
for cell = 1:nc
    incell = find(zK >= z(cell) & zK<z(cell+1));
    if ~isempty(incell)
        Kc=K(incell);w=ones(size(Kc));
        Ks(cell) = vavg(Kc,w,avgmeth);
    else
        %         disp(strcat(['...>>> cell: ',num2str(cell),' centered at ',...
        %                      num2str(zm(cell)),' m empty']))
    end
end

uppval=find(isfinite(Ks),1,'first');
if ~exist('KTop','var'), KTop=Ks(uppval);end
Ks(1:uppval-1)=KTop;

lowval=find(isfinite(Ks),1,'last');
if ~exist('KBot','var'),KBot=Ks(lowval);end
Ks(lowval+1:nc)=KBot;

okval=find(isfinite(Ks));zi=zm(okval);Ki=Ks(okval);
Ks=interp1(zi,Ki,zm);


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
    plot(K,zK,'.b','LineWidth',1); hold on
    plot(Ks,zm,'-r','LineWidth',2); hold on
    
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
    Kx=Ks(id);Kx=Kx(1:length(Kx)-1)';
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

k=Ks;
kA=0.0013*nones';
kB=0.0029*nones'; % OKU fit Kukkonen
%kA=0.00*nones';
%kB=0.00*nones';


r=2800*nones';
c=740*nones';
%rc=r.*c;
rc=rcl(k);  % OKU fit Kukkonen

h = 1.0e-6*nones';
p = 0.01*nones';

z_layer=[  0,  33,  1314,  1515,  1658,  1677,  1863, 1900, 2001, 2213,  2304, 2415, 3000, 5000];
h_sedim = 1.00e-6;
h_metam = 1.70e-6;
h_ophio = 1.55e-6;
h_pegma = 6.*1.e-6;
h_layer=[h_sedim   h_metam  h_ophio,  h_metam    h_pegma   h_metam   h_pegma   h_metam   h_pegma   h_metam  h_pegma  h_metam h_pegma];
p_cryst = 0.02;
p_sedim = 0.2;
p_layer=[  p_sedim p_cryst p_cryst  p_cryst  p_cryst  p_cryst  p_cryst p_cryst  p_cryst  p_cryst p_cryst  p_cryst p_cryst];


for ih=1:length(z_layer)-1
    layer=find(z>=z_layer(ih) & z<z_layer(ih+1));
    h(layer) = h_layer(ih);
    p(layer) = p_layer(ih);
end

h(zm>zDatBot)=0.;

% STRUCTURE DATA
nd=length(id);
errT=ErrDeflt;
Terr=errT*ones(size(id))';
Tcov=spdiags(Terr.^2,0,nd,nd);


sitepar=mstruct(k,kA,kB,h,p,c,r,rc,z,ip,t,it,qb,gts,Tobs,id,zobs,Tcov,Terr,props,name);


F=strcat([name '_SitePar.mat']);
save(F,'sitepar');
disp([' >>>>> site parameter saved to:' F]);
disp([' ']);


end
