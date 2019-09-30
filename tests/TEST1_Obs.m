function [ierr]=SITE_Obs(name)
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

Data_file 		= 'DataBallingB.csv';
site            = 'TEST1';
props           = 'const';
out             = 0;
estq            = 0;
estq_int        = [100 200 300 400 500];
plotit          = 0;

% PARAMETER FOR PREPROCESSING
zDatTop         =   10;
zDatBot         =   2000;

GST0            =   0;
POM             =   0; 
POR             =   0.00;
KTop            =   2.5e+00;   % mean aus TC-m top 100
KBot            =   2.5e+00; % mean aus TC-m entire well
KMean           =   2.5e+00;
ErrDeflt        =   0.1;

L=3;
CovType = 'g';

% s='shuffle';
% s=rng;

smooth_props='m';
avgmeth ='h';
smooth_data='s';
nspline=101;
wspline=0.5;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE sITE_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_Obs_in.mat']);
if exist(F)
    disp([mfilename ' defaults overwritten from ', F])
    load(F); mstruct(prep_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if plotit
    set_graphpars
end

disp(['   ']);
disp(strcat([ ' ...Preprocessing site ', site]));

step = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: GET SPATIAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': set grids']));
% SPATIAL AND TEMPORAL MESH
disp([' ...set grids'])

meshfileZ=[name,'_DepthGrid.mat'];
load (meshfileZ);

% disp([mfilename '    spatial mesh: ',num2str(nz),' temporal mesh:',num2str(nt)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: READ & PREPROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': read obs']));

% 
% O1    =   importdata([datpath,Data_file]);
% 
% % Version 2012b
% zT = O1(:,1); 
% T = O1(:,4); zT = zT(isfinite(T)); T = T(isfinite(T));
% 
% N = length(T(:,1));
% 
% if L~= 0
%     if exist('s')
%         rng(s);
%     end
%     
%     switch lower(CovType)
%         case{'e' 'markov' 'exponential'}
%             Cov=CovarExpnl(ones(N,1),L);
%         case{'g' 'gauss'}
%             Cov=CovarGauss(ones(N,1),L);
%             %     case{'m' 'matern'}
%             %         Cov=CovarMatern(ones(N,1),L);
%         otherwise
%             error([lower(CovFun), ' not implemented. STOP']);
%     end
%     C    = chol(Cov);
%     
%     Err = ErrDeflt;
%     
%     err_nor =  Err*randn(N,1);
%     err_cor =  err_nor'*C;
%     
%     T = T+err_cor';
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% step=step+1;
% disp(strcat([ ' ...>>> Step ',num2str(step),...
%     ': average & interpolate to FD mesh']));
% 
% 
% % CELL CENTERS
% zm=0.5*(z(1:nz-1)+z(2:nz));nc=length(zm);
% 
% % TEMPERATURE
% breaks=linspace(min(zT),max(zT),nspline);
% ss=splinefit(zT,T,breaks,'r',wspline);dd=ppdiff(ss);
% Tx=ppval(ss,z);Tc=ppval(ss,zm);dTc=ppval(dd,zm);
% resT=ppval(ss,zT)-T;no=length(resT);errT=norm(resT)/sqrt(no);
% Ts=NaN(size(z));
% Ts(z>=min(zT) & z <=max(zT))=Tx(z>=min(zT) & z <=max(zT));
% 
% 
% % CUT LOG
% Ts(z< zDatTop)=NaN;
% Ts(z> zDatBot)=NaN;
% id=find(isfinite(Ts));
% Tobs=Ts(id);zobs=z(id); nd=length(id);Tobs=Tobs';
% 
% if plotit
%     % PLOT TEMPERATURES
%     figure
%     plot(T,zT,'.b','LineWidth',1); hold on
%     plot(Tobs,zobs,'-r','LineWidth',1); hold on
% %     xlim([0 50]);
% %     ylim([0 2750]);
%     set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
%     xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
%     ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
%     title(strcat([site, ' temperatures']),'FontSize',fontsz,'FontWeight',fontwg);
%     S1=strcat(['orig']);
%     S2=strcat(['mod']);
%     legend(S1,S2,'location', 'southwest')
%     grid on
%     file=strcat([name '_Temp']);
%     saveas(gcf,file,plotfmt)
%     
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP IV: DEFINE STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ ' ...>>> Step ',num2str(step),': setup structure']));
%

% STRUCTURE DATA
% nd=length(id);
% errT=ErrDeflt;
% Terr=errT*ones(size(id))';
% Tcov=spdiags(Terr.^2,0,nd,nd);
id          = []; nd=length(id);
zobs        = [];
Tobs        = [];
Terr        = [];    
Tcov        = [];  



siteobs=mstruct(Tobs,id,zobs,Tcov,Terr,name);


F=strcat([name '_SiteObs.mat']);
save(F,'siteobs');
disp([' >>>>> site observations saved to:' F]);
disp([' ']);


end
