% PLOT OF INVERSIONS
% clear all
% close all
% clc

load common

debug=0;

% SET PATHS
addpath([srcpath,'/tools/'])
addpath([srcpath,'/src/'])

dfmt=1;ffmt='.zip';
archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);
SITE='TMP';

% LOCAL OPTIONS
yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS
plotfmt1='png';
plotfmt2='epsc2';
outsteps=false;
graph=false;
fontwg='normal';
fontsz=16;
cols={'-b','-r','-g','-m','-k' '--b','--r','--g','--m','--k' ':b',':r',':g',':m',':k'};

plot_resd=1;
plot_gsth=1;
plot_grad=1;
plot_tlog=1;
include_proxies=0;
smooth_props='m';
avgmeth ='h';
smoothpl=1;

smooth_data='m';
nsplineD=101;wsplineD=0.5;
smoothdl=1;

smooth_grads = 'm';
smoothgl=1;
gmode={'box','mir'};

% 
% if include_proxies
%     load([ datpath 'Luterbacher2004_FMK']);
%     tLut=abs(YL-2002)+2;
%     TLut=0.71*T30L+2.9; %ILmo
%     TLut=0.59*T30L+2.76; %Jan
%     lLut=strcat(['Luterbacher 2004']);
% end

Tlimits=[-15 20];

zlimits=[0 4000];
if plot_tlog,
    tscal=1.;
    %tlimits=[10 150000];
    tlimits=[10 150000];
    
else
    tscal=1e-3;
    tlimits=[0 16];
end

disp(['   '])
run='GCV - Q_b Experiment';
%run = strrep()
TMP=dir('*_results.mat');
Files={TMP.name};
%Files={...
% };
   


Nfiles=length(Files);
whichfiles=[1:Nfiles];
nfiles=length(whichfiles);

INFO=fopen('INFO.dat','w');
imodel=0;
for ifiles=1:nfiles
    imodel=imodel+1;
    filename=Files{ifiles};
    %disp(['   ']);disp([filename,' rms =',num2str]);
    load(filename)
    %reg1par
    allt{imodel}=t;allit{imodel}=it;
    [nit,nst]=size(m_iter);
    allname{imodel}=name;
    mm=m_iter(nit,:);mm=mm(it);
    allres{imodel}=r_iter(nit,:)';
    allgsth{imodel}=mm(:)';
    allrms{imodel}=rms_iter(nit);
    allq0{imodel}=qb;
    allt0{imodel}=gts;
    allmeth{imodel,1}=0;
    allmeth{imodel,2}=0;
    allmeth{imodel,3}=regpar0;
    allcov{imodel}=Cmm;

    smooth_grads='m';smoothgl=21;gmode={'box','mir'};M=1;
    zd=z(id);dz=diff(z(id));ii=id(1:length(id)-1)';
    zm=0.5*(zd(1:length(zd)-1)+zd(2:length(zd)));
    TC=Tcalc(id,nt);dT=diff(TC);dTC=dT(:)./dz(:);QC = Tcon(ii).*dTC(:);
    TO=Tobs;dT=diff(TO);dTO=dT(:)./dz(:);Q0 = Tcon(ii).*dTO(:);
    switch lower(smooth_grads)
        case {'s','spline'}
            spoints=linspace(min(zG),max(zG),nsplineG);
            pp=Q0(:);ss=splinefit(zG,pp,spoints,'r',wsplineG);
            Qobs=ppval(ss,zG);
        case {'f','filter','m','mavg'}
            pp=Q0(:);npp=length(pp);LPROP=smoothgl;
             ss= wfilt(pp,LPROP,M,gmode);
            Qobs=ss;
        otherwise
            disp([' no smoothing!'])
    end
    
    
    allTcal{imodel}=TC(:)';
    allTobs{imodel}=TO(:)';
    allQcal{imodel}=QC(:)';
    allQobs{imodel}=Qobs(:)';
    allTCon{imodel}=Tcon(:)';
    allzm{imodel}=zm(:)';
    allzd{imodel}=zd(:)';
    disp(['   ']);disp([filename,',  rms =',num2str(allrms{imodel})]);

end

%
% PLOT RESULTS
disp(['   ']);disp([' ...plot results']);
%
% FIG GSTH

[n1,n2]=size(allgsth);
N=length(allres);
%

icurv=0;
for ii=1:N
    icurv=icurv+1;
%     %regpval=allmeth{ii,3};sregpval=num2str(regpval(1));
    rmsval = allrms{ii};
    q0val  = allq0{ii};
    legstr{icurv}= ...
         strcat([num2str(abs(q0val*1000)),'/',num2str(rmsval,'%5.3f')]);

%     %     legstr{icurv}= ...
%     %         strcat([sregpval,'/',num2str(allrms(ii),'%5.3f'),'/',allmeth{ii,2}]);
%     
%     %     legstr{icurv}= ...
%     %         strcat([run ': ',num2str(allq0(ii)*1000),'/',num2str(allrms(icurv),'%3.2f')]);
%     legstr{icurv}= ...
%         strcat([run ': ',QN{ii},'/',num2str(allrms{icurv},'%3.2f')]);
%     
end
% GSTH
if plot_gsth
    
    
    figure;
    
    if plot_tlog
        for ii=1:N
               ty=allt{ii}/year2sec+13.5;
            h(ii)=plot(-ty*tscal,allgsth{ii},cols{ii},'LineWidth',2);hold on
        end
%         if include_proxies
%             legstr=[legstr lLut];
%             plot(tLut,TLut,'+b','LineWidth',1);
%             %legstr=[legstr lLuto];
%             %plot(tLut,TLut_orig,'.m','LineWidth',1);
%         end
        %xlim(tlimits);
        set(gca,'xscale','log','xdir','rev',...
            'xtick',[10 100 1000 10000 100000],...            
            'FontSize',fontsz,'FontWeight',fontwg);
        locl='southeast';
    else
        for ii=1:N
            h(ii)=plot(-ty*tscal,allgsth{ii},cols{ii},'LineWidth',2);hold on
        end
%         if include_proxies
%             legstr=[legstr lLut];
%             plot(tLut,TLut,'+b','LineWidth',1);
%             %plot(tHist,THist,'*m','LineWidth',1);
%         end

        %xlim(tlimits);
        set(gca,'xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        %            'xtick',[0:2:16],...
        locl='southeast';
    end
    title(strcat([SITE,': ',run]),'FontSize',fontsz,'FontWeight',fontwg);
    legend(legstr,'location','best');
    grid on;
    xlabel('Time BP (a)');ylabel('T (C)');
    filename=strcat(SITE,'_GSTH_',run);
    saveas(gcf,filename,plotfmt1)
        saveas(gcf,filename,plotfmt2)
end


%RESIDUALS
if plot_resd
    figure;
    for ii=1:N    
        plot(allres{ii},allzd{ii},cols{ii},'LineWidth',2);hold on
    end
    grid on;
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('depth (m)');xlabel('\Delta T (K)');
    %ylim(zlimits);
%     %xlim([-0.225 0.075]);
    %%xlim([-0.02 0.02]);
    title(strcat([SITE,': ',run]),'FontSize',fontsz,'FontWeight',fontwg);
    legend(legstr,'location','northwest');
    filename=strcat(SITE,'_RESID_',run);
    saveas(gcf,filename,plotfmt1)
        saveas(gcf,filename,plotfmt2)
end

if plot_grad,
    % GRADIENTS
    figure
    plot(allQobs{1}*1000,allzm{1},'-k','LineWidth',3); hold on
    for ii=1:N
        plot(allQcal{ii}*1000,allzm{ii},cols{ii},'LineWidth',2); hold on
    end
     %ylim(zlimits);
%     %xlim([20 60]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('\lambda \nabla T (K/km)','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat([SITE,': ',run]),'FontSize',fontsz,'FontWeight',fontwg);
    obsstr = {'observed'};
    str = [obsstr legstr];
    legend(str,'location','southwest')
    grid on
    filename=strcat([SITE '_Q0',]);
    %export_fig filename -transparent -png   %-pdf -eps
    saveas(gcf,filename,plotfmt1)
        saveas(gcf,filename,plotfmt2)
end

