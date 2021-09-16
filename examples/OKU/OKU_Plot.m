
% PARAMETERS FOR INVERSION
clear all
close all
clc
warning('off','all')

yr2s=31557600;s2yr=1/yr2s;


load('common.mat')

addpath(['local'])
addpath([srcpath,'./src']);
addpath([srcpath,'./src/mcmc']);
addpath([srcpath,'./tools'])


disp(['   '])

nburnin=1000;
Terr=0.1;Terr2=Terr^2;

site  ='OKU';props  =lower(site);meth='DR';verstrng='_Z';
addpath([srcpath,strcat(['src/props/',lower(props)])])

disp([' ']);disp(strcat([' plot MCMC results for ' site ]));


TMP=dir('*save.mat');
SaveFiles={TMP.name};

% TMP1=dir('LAX_TikhFIX*results.mat');
% Files1={TMP1.name};
% TMP2=dir('LAX_TikhU*results.mat');
% Files2={TMP2.name};
% Files=[Files1 Files2];





%SaveFiles=get_contents(pwd,'filter','save.mat');
%SaveFiles={'OLK_MCMC_15-Nov-2014_92208_save.mat'};
Nfiles=length(SaveFiles);
whichfiles=[1:Nfiles];
%whichfiles=[ 2 3 4  6];
%whichfiles=[ 1 2 3 ];

nfiles=length(whichfiles);

GSTH=[];QB=[];HB=[];RZ=[];RMS=[];TC=[];nrej=0;

for ifiles=1:nfiles
    File=SaveFiles{whichfiles(ifiles)};
    disp([' ']);disp(strcat([' read MCMC chain #',num2str(ifiles),' from ' File ]));
    filename=strcat([File]);
    load(filename)
    [n1 n2]=size(parchain);
    GSTHX=parchain(nburnin:n1,1:n2-2);
    FILL=isfinite(GSTHX(:,1));GSTHX=GSTHX(FILL,1:n2-2);
    GSTH=[GSTH;GSTHX];
    QBX=parchain(nburnin:n1,n2-1);QBX=QBX(FILL);
    QB=[QB;QBX];
    HX=parchain(nburnin:n1,n2-0);HX=HX(FILL);
    HB=[HB;HX];
    RX=reschain(nburnin:n1,:);    RX=RX(FILL,:);
    RZ=[RZ;RX];
    CX=calchain(nburnin:n1,:);    CX=CX(FILL,:);
    TCZ=[TC;CX];
    [n1 n2]=size(RX);
    RMSX=sqrt(sum(RX.^2./Terr2,2)/((n2-1)));
    RMS=[RMS;RMSX];
    nrej=nrej+rej;
end
% whos QB GSTH RMS




MC=GSTH;


[nsmp,npars]=size(GSTH);
nsimu=nsmp;

sitepar=data.sitepar;mstruct(sitepar);
numpar =data.numpar;mstruct(numpar);


S=site; N=site;
%
% disp([' ']);disp(strcat([' generate model for ' N ]));
% CMD=strcat(['Prep_' S  verstrng '(S,N);']);eval(CMD);


% disp([' ']);disp([' >>>>> site configuration read']);
Z=z(id);Z=Z(:)';

Zedges=[0 0.5*(Z(1:length(Z)-1)+Z(2:length(Z)))  Z(length(Z))+Z(length(Z))-Z(length(Z)-1)];
ty=-t*s2yr;


%load(strcat(['OLKMCMC_prior']));
nsteps=25;
m_apr=[5*ones(nsteps,1)];
GST_PRIOR=m_apr(it);nsteps=length(m_apr);


Q=[.68  .95];P1=[1-Q(1) Q(1)];P2=[1-Q(2) Q(2)];
S=[  1    2];


Rmed=median(RMS);Ravg=mean(RMS);Rstd=std(RMS);Rmad=mad(RMS);
[R1]=quantile(RMS,P1,1);
[R2]=quantile(RMS,P2,1);
disp([' RMS mean = ',num2str(Ravg) ', std = ',num2str(Rstd)]);
disp([' RMS med  = ',num2str(Rmed) ', Q1/2 = ',num2str(R1(1)),'/',num2str(R1(2))]);

Qmed=median(QB);Qavg=mean(QB);Qstd=std(QB);Qmad=mad(QB);
[Q1]=quantile(QB,P1,1);
[Q2]=quantile(QB,P2,1);
disp([' HFD mean = ',num2str(Qavg),',  std = ',num2str(Qstd)]);
disp([' HFDmed  = ',num2str(Qmed),', Q1/2 = ',num2str(Q1(1)),'/',num2str(Q1(2))]);

Hmed=median(HB);Havg=mean(HB);Hstd=std(HB);Hmad=mad(HB);
[R1]=quantile(HB,P1,1);
[R2]=quantile(HB,P2,1);
disp([' HProd mean = ',num2str(Havg) ', std = ',num2str(Hstd)]);
disp([' HProd med  = ',num2str(Hmed) ', Q1/2 = ',num2str(R1(1)),'/',num2str(R1(2))]);


%GRAPHICS
CovFun='G';L=5;
NAME=strcat([site, '_',meth,CovFun,num2str(L),'_GQ']);
tit = NAME;



Tlimits=[-4 10];
tlimits=[10 150000];
zlimits=[0  1500];
rlimits=[-0.1  0.1];

plot_reslt=1;
plot_qhist=1;
plot_hhist=1;
plot_rhist=1;
plot_tresd=1;
plot_gresd=1;
plot_tgrad=1;
plot_med = 1;
plot_avg= 0;
plot_cmap=0;

include_proxies=0;
include_levmodl=0;
include_prior=1;
include_text=1;

set_graphpars

% plotfmt='epsc2'
%==========================================================================
if plot_reslt
    
    TY=abs(t*s2yr)-13.5+eps;
    
    if plot_med
        
        % Median and Quantiles
        GC=nanmedian(MC,1);
        [GQ1]=quantile(MC,P1,1);
        [GQ2]=quantile(MC,P2,1);
        GC=GC(it); GQ1=GQ1(:,it);GQ2=GQ2(:,it);
        legstr={[num2str(Q(2)*100),'% quantiles'],...
            [num2str(Q(1)*100),'% quantiles'],'median'};
        
    else
        % Mean and std
        GC=nanmean(MC,1); GS=nanstd(MC,1);
        GQ1=[GC(:)-S(1)*GS(:);GC(:)+S(1)*GS(:)]';
        GQ2=[GC(:)-S(2)*GS(:);GC(:)+S(2)*GS(:)]';
        GC=GC(it); GQ1=GQ1(:,it);GQ2=GQ2(:,it);
        legstr={[num2str(S(1)),'*sd'],[num2str(S(2)),'*sd'],'average'};
    end
    %--------------------------------------------------------------------------
    save('LX2_median','TY','GC','GQ1','GQ2'),
    
    figure;
    
    
    if shaded
        cropsh=TY > tlimits(1)+10*eps & TY< tlimits(2)-10*eps;
        TX=TY(cropsh);GQ1=GQ1(:,cropsh);GQ2=GQ2(:,cropsh);
        
        %plot(TY,flipud(GC),'-r','LineWidth',3);hold on
        TX=TX(:);TX=[TX;flipud(TX)];NX=length(TX);
        GQa=GQ2(1,:); GQa=GQa(:);
        GQb=GQ2(2,:); GQb=GQb(:) ;
        GQX=[GQa;flipud(GQb)];
        fill(TX,GQX,[grey1 grey1 grey1],'edgecolor',shedge);hold on
        GQa=GQ1(1,:); GQa=GQa(:);
        GQb=GQ1(2,:); GQb=GQb(:) ;
        GQX=[GQa;flipud(GQb)];
        fill(TX,GQX,[grey2 grey2 grey2],'edgecolor',shedge);hold on
        plot(TY,GC,'-r','LineWidth',3);hold on
        
    else
        plot(TY,GC,'-r','LineWidth',3);hold on
        plot(TY,GQ1,'--b','LineWidth',1);hold on
        plot(TY,GQ2,':b','LineWidth',1);hold on
        
    end
    
    grid on;
    set(gca,'xscale','log','xdir','rev',...
        'xtick',[10 100 1000 10000 100000],...
        'ytick',[-4 0 4 8],...
        'FontSize',fontsz,'FontWeight',fontwg);
    xlabel('Time BP/2000 (a)');ylabel('\Delta T (K)');
    ylim(Tlimits);
    xlim(tlimits);
    %title([site,': GSTH, MCMC average'],'FontSize',fontsz)
    %legend ('MCMC average','MCMC Median','MCMC Prior','location','southeast')
    filename=[tit '_GSTH1_N',num2str(nsmp,'%i')];
    if include_prior
        filename=[filename '_prior'];
        legstr=[legstr 'prior'];
        plot(TY,GST_PRIOR,'-k','LineWidth',3); hold on
    end
    if include_proxies
        load([ datpath 'Luterbacher2004_LAX']);
        tLut=abs(YL-2000);
        TLut0=T30L;
        TLut1=0.71*T30L+2.93; % FINLAND
        TLut2=0.59*T30L+2.76; % FORSMARK
        filename=[filename '_proxies'];
        legstr=[legstr 'Luterbacher 2004' 'Heikkilae 2003'];
        %        plot(tLut,TLut0,'+b','LineWidth',1);
        %         plot(tLut,TLut1,'dm','LineWidth',1);
        plot(tLut,TLut1,'xb','LineWidth',1);
        %plot(tHei,THei,'*m','LineWidth',1);
    end
    if include_levmodl
        filename=[filename '_glcmod'];
        legstr=[legstr 'GSM samples'];
        plot(tLev,TLev,':','LineWidth',2);
    end
    if include_text
        dtit=strrep(tit,'_',' ');
        dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
            num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
        textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    end
    legend(legstr,'location','southeast');
    saveas(gcf,filename,plotfmt)
end


if plot_qhist
    Qmed=median(QB);Qavg=mean(QB);Qstd=std(QB);Qmad=mad(QB);
    [Q1]=quantile(QB,P1,1);
    [Q2]=quantile(QB,P2,1);
    disp('HFD');
    disp([' mean = ',num2str(Qavg),',  std = ',num2str(Qstd)]);
    disp([' med  = ',num2str(Qmed),', Q1/2 = ',num2str(Q1(1)),'/',num2str(Q1(2))]);
    
    %--------------------------------------------------------------------------
    figure
    %     QB(QB< 20) =NaN;
    %     QB=QB(10000:length(QB));
    hist(QB,32);
    grid on
    %title([site,': Q_{base} posterior, MCMC '],'FontSize',fontsz);
    xlabel('HFD (mW/m^2)','FontSize',fontsz);ylabel('N (-)','FontSize',fontsz);
    %xlim([37 43]);
    if include_text
        dtit=strrep(tit,'_',' ');
        dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
            num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
        textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    end
    set(gca,'FontSize',fontsz,'FontWeight',fontwg);
    
    filename=strcat(tit,'_Q0');
    saveas(gcf,filename,plotfmt)
    %--------------------------------------------------------------------------
end

if plot_rhist
    
    Rmed=median(RMS);Ravg=mean(RMS);Rstd=std(RMS);Rmad=mad(RMS);
    [R1]=quantile(HB,P1,1);
    [R2]=quantile(HB,P2,1);
    disp('RMS');
    disp([' mean = ',num2str(Ravg) ', std = ',num2str(Rstd)]);
    disp([' med  = ',num2str(Rmed) ', Q1/2 = ',num2str(R1(1)),'/',num2str(R1(2))])
    figure
    hist(HB,32);
    grid on
    %title([site,': Q_{base} posterior, MCMC '],'FontSize',fontsz);
    xlabel('H (\mu W/m^{-3})','FontSize',fontsz);ylabel('N (-)','FontSize',fontsz);
    
    if include_text
        dtit=strrep(tit,'_',' ');
        dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
            num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
        textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    end
    set(gca,'FontSize',fontsz,'FontWeight',fontwg);
    filename=strcat(tit,'_HProd');
    saveas(gcf,filename,plotfmt)
    %--------------------------------------------------------------------------
end

if plot_hhist
    
    Hmed=median(HB);Havg=mean(HB);Hstd=std(HB);Hmad=mad(HB);
    [H1]=quantile(HB,P1,1);
    [H2]=quantile(HB,P2,1);
    disp('H');
    disp([' mean = ',num2str(Havg) ', std = ',num2str(Hstd)]);
    disp([' med  = ',num2str(Hmed) ', Q1/2 = ',num2str(H1(1)),'/',num2str(H1(2))]);
    
    figure
    hist(RMS,32);
    grid on
    %title([site,': Q_{base} posterior, MCMC '],'FontSize',fontsz);
    xlabel('RMS (-)','FontSize',fontsz);ylabel('N (-)','FontSize',fontsz);
    
    if include_text
        dtit=strrep(tit,'_',' ');
        dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
            num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
        textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    end
    set(gca,'FontSize',fontsz,'FontWeight',fontwg);
    filename=strcat(tit,'_RMS');
    saveas(gcf,filename,plotfmt)
    %--------------------------------------------------------------------------
end

%==========================================================================
if plot_tresd
    
    
    if plot_med
        % Median and Quantiles
        RC=nanmedian(RZ,1);
        [RQ1]=quantile(RZ,P1,1);
        [RQ2]=quantile(RZ,P2,1);
        legstr={'median',[num2str(Q(2)*100),'% quantiles'],...
            [num2str(Q(1)*100),'% quantiles']};
        
    else
        % Mean and std
        RC=nanmean(RZ,1); Tstd=nanstd(RZ,1);
        RQ1=[RC(:)-S(1)*Tstd(:);RC(:)+S(1)*Tstd(:)]';
        RQ2=[RC(:)-S(2)*Tstd(:);RC(:)+S(2)*Tstd(:)]';
        legstr={'average', [num2str(S(1)),'*sd'],[num2str(S(2)),'*sd']};
    end
    %
    
    
    figure;
    
    if shaded
        plot(RC,Z/1000,'-r','LineWidth',3);hold on
        ZX=Z/1000;ZX=ZX(:);ZX=[ZX;flipud(ZX)];
        
        RQa=RQ2(1,:); RQa=RQa(:);
        RQb=RQ2(2,:); RQb=RQb(:) ;
        RQX=[RQa;flipud(RQb)];
        fill(RQX,ZX,[grey1 grey1 grey1]);
        RQa=RQ1(1,:); RQa=RQa(:);
        RQb=RQ1(2,:); RQb=RQb(:) ;
        RQX=[RQa;flipud(RQb)];
        fill(RQX,ZX,[grey2 grey2 grey2]);
        plot(RC,Z/1000,'-r','LineWidth',2);hold on
        legend(legstr,'location','west');
        
    else
        plot(RC,Z/1000,'-r','LineWidth',2);hold on
        plot(RQ1,Z/1000,'--b','LineWidth',1);hold on
        plot(RQ2,Z/1000,':b','LineWidth',1);hold on
        legend(legstr,'location','west');
        
    end
    grid on
    xlim(rlimits);
    ylim(zlimits/1000);
    ylabel('z (km)','FontSize',fontsz);
    
    xlabel('T (C)','FontSize',fontsz);
    set(gca,'YDir','reverse')
    set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
    dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
        num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
    textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    
    file=strcat([tit '_Residuals']);
    saveas(gcf,file,plotfmt);
    
    
end


%==========================================================================
if plot_gresd
    [n1,n2]=size(RX);
    ZM=0.5*(Z(1:n2-1)+Z(2:n2));
    GX=diff(RX,1,2)./repmat(diff(Z),n1,1);
    KX=mean(k)*ones(size(GX));
    QX=KX.*GX;
    ZM=0.5*(Z(1:length(Z)-1)+Z(2:length(Z)));
    depthint=ZM>100 & ZM <1350;
    ZM=ZM(depthint); GX=GX(:,depthint);
    
    
    if plot_med
        % Median and Quantiles
        GC=nanmedian(GX,1);
        [GQ1]=quantile(GX,P1,1);
        [GQ2]=quantile(GX,P2,1);
        legstr={'median',[num2str(Q(2)*100),'% quantiles'],...
            [num2str(Q(1)*100),'% quantiles']};
        
    else
        % Mean and std
        GC=nanmean(GX,1); Tstd=nanstd(GZ,1);
        GQ1=[GC(:)-S(1)*Tstd(:);GC(:)+S(1)*Tstd(:)]';
        GQ2=[GC(:)-S(2)*Tstd(:);GC(:)+S(2)*Tstd(:)]';
        legstr={'average', [num2str(S(1)),'*sd'],[num2str(S(2)),'*sd']};
    end
    
    
    
    figure;
    
    if shaded
        
        plot(GC,ZM/1000,'-r','LineWidth',3);hold on
        ZM=ZM(:)/1000;ZX=[ZM;flipud(ZM)];
        
        GQa=GQ2(1,:); GQa=GQa(:);
        GQb=GQ2(2,:); GQb=GQb(:) ;
        GQX=[GQa;flipud(GQb)];
        fill(GQX,ZX,[grey1 grey1 grey1]);
        GQa=GQ1(1,:); GQa=GQa(:);
        GQb=GQ1(2,:); GQb=GQb(:) ;
        GQX=[GQa;flipud(GQb)];
        fill(GQX,ZX,[grey2 grey2 grey2]);
        plot(GC,ZM,'-r','LineWidth',2);hold on
        legend(legstr,'location','west');
        
    else
        plot(GC,ZM,'-r','LineWidth',2);hold on
        plot(GQ1,ZM,'--b','LineWidth',1);hold on
        plot(GQ2,ZM,':b','LineWidth',1);hold on
        legend(legstr,'location','west');
        
    end
    grid on
    %xlim(rlimits);
    ylim(zlimits/1000);
    ylabel('z (km)','FontSize',fontsz);
    xlabel('\delta T (K)','FontSize',fontsz);
    
    set(gca,'YDir','reverse')
    set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
    dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
        num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
    textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    file=strcat([tit '_GResiduals']);
    saveas(gcf,file,plotfmt);
    
    
end

%==========================================================================
if plot_tgrad
    shaded=0;
    [n1,n2]=size(CX);
    
    ZM=0.5*(Z(1:n2-1)+Z(2:n2));
    TO=sitepar.Tobs;
    GO=diff(TO,1)'./diff(Z,1);
    
    LF=11;TF={'tri','mir'};M=1;
    TO= wfilt(TO,LF,M,TF);
    GO=diff(TO,1)'./diff(Z,1);
    
    
    
    GX=diff(CX,1,2)./repmat(diff(Z),n1,1);
    
    %KX=mean(k)*ones(size(GX));
    %QX=KX.*GX;
    
    depthint=ZM>100 & ZM <1400;
    ZM=ZM(depthint); GX=GX(:,depthint); GO=GO(:,depthint);
    
    legstr=[];
    legstr=[legstr,{'observed'}];
    if plot_med
        % Median and Quantiles
        GC=nanmedian(GX,1);
        [GQ1]=quantile(GX,P1,1);
        [GQ2]=quantile(GX,P2,1);
        legstr=[legstr {'median'},{[num2str(Q(2)*100),'% quantiles']},...
            {[num2str(Q(1)*100),'% quantiles']}];
    else
        % Mean and std
        GC=nanmean(GX,1); Tstd=nanstd(GZ,1);
        GQ1=[GC(:)-S(1)*Tstd(:);GC(:)+S(1)*Tstd(:)]';
        GQ2=[GC(:)-S(2)*Tstd(:);GC(:)+S(2)*Tstd(:)]';
        legstr=[legstr, {'average'}, {[num2str(S(1)),'*sd']},{[num2str(S(2)),'*sd']}];
    end
    
    
    
    figure;
    ZM=ZM(:)/1000;
    
    plot(GO,ZM,'-+b','LineWidth',1);hold on
    plot(GC,ZM,'-r','LineWidth',2);hold on
    plot(GQ1,ZM,'--r','LineWidth',1);hold on
    plot(GQ2,ZM,':r','LineWidth',1);hold on
    
    legend(legstr,'location','southwest');
    
    grid on
    %xlim(rlimits);
    ylim(zlimits/1000);
    ylabel('z (km)','FontSize',fontsz);
    xlabel('grad T (K/m)','FontSize',fontsz);
    set(gca,'YDir','reverse')
    set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
    dtext= strcat([dtit,': ',num2str(nsmp,'%i'),' samples, ',...
        num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
    textloc(dtext,'north','FontSize',fontsz-4,'FontWeight',fontwg);
    
    file=strcat([tit '_GradientT']);
    
    saveas(gcf,file,plotfmt);
    
    
end


%
% %==========================================================================
if plot_cmap
    %
    %     ty=t*s2yr;
    %     tstart=tstart*s2yr;
    %     tend  =tend *s2yr;
    %     tedges = logspace(log10(tstart),log10(tend),nsteps); tedges=[tedges 0.1*tedges(nsteps)];
    %     tedges=fliplr(tedges);
    %     tc=exp(0.5*(log(tedges(1:ynsteps))+log(tedges(2:nsteps+1))));
    %     Tedges = linspace(-10,14,120);;
    %     [n1 n2]=size(GSTH);GSTHV=fliplr(GSTH);
    %     GSTHV=reshape(GSTHV',n1*n2,1);%GSTHV=fliplr(GSTHV);
    %     TIMEV= repmat(tc(:),n1,1);
    %     histmat = histnd(TIMEV,GSTHV, tedges, Tedges);
    %     histmat=histmat/max(max(histmat));
    %
    %     %--------------------------------------------------------------------------
    %     figure;
    %     pcolor(tedges,Tedges,histmat'); colorbar ;hold on
    %     %shading interp
    %     colormap(colscl);
    %     shading(gca,shade)
    %     xlim(tlimits);
    %     ylim(Tlimits);
    %     set(gca,'xscale','log','xdir','rev','FontSize',fontsz,'FontWeight',fontwg);
    %     set(gca,'tickdir','out','FontSize',fontsz,'FontWeight',fontwg);
    %     %title([site,': GSTH, MCMC ',num2str(nsmp,'%.1e'),' samples'],'FontSize',fontsz)
    %     xlabel('Time BP (a)');ylabel('\Delta T (K)');
    %     filename=strcat(tit,'_GSTH1');
    %     saveas(gcf,filename,plotfmt)
    %     plot(abs(ty),gstp,'-k','LineWidth',2)
    %     plot(tLev,TLev,':w');
    %     filename=strcat(tit,'_GSTH1prior');
    %     sav    disp('H');eas(gcf,filename,plotfmt)
    %
    %     %--------------------------------------------------------------------------
    %     figure;
    %     pcolor(tedges,Tedges,histmat'); colorbar ;hold on
    %     colormap(colscl);
    %     shading(gca,shade)
    %     xlim(tlimits);
    %     ylim(Tlimits);
    %     set(gca,'xscale','log','xdir','rev','FontSize',fontsz,'FontWeight',fontwg);
    %     set(gca,'tickdir','out','FontSize',fontsz,'FontWeight',fontwg);
    %     %title([site,': GSTH, MCMC ',num2str(nsmp,'%.1e'),' samples'],'FontSize',fontsz)
    %     xlabel('Time BP (a)');ylabel('\Delta T (K)');
    %     plot(tLut,TLut,'+b','LineWidth',1);
    %     plot(tHei,THei,'*m','LineWidth',1);
    %     plot(tLev,TLev,':w');
    %     plot(abs(ty),gstp,'-k','LineWidth',2);
    %     filename=strcat(tit,'_GSTHp1roxy');
    %     saveas(gcf,filename,plotfmt)
    %
    %     %--------------------------------------------------------------------------
    %     figure;
    %     plot(abs(ty),gstp,'-k','LineWidth',3); hold on
    %     plot(abs(ty),gsta,'-r','LineWidth',3);hold on
    %     plot(tLut,TLut,'+b','LineWidth',1)
    %     plot(tHei,THei,'*m','LineWidth',1)
    %     plot(tLev,TLev,':','LineWidth',2)
    %     grid on
    %     xlim(tlimits);
    %     ylim(Tlimits);
    %     set(gca,'xscale','log','xdir','rev','FontSize',fontsz,'FontWeight',fontwg);
    %     set(gca,'tickdir','out','FontSize',fontsz,'FontWeight',fontwg);
    %     dtext= strcat([site,': MCMC ',num2str(nsmp,'%i'),' samples, ',...
    %         num2str((nsmp-nrej)*100/nsmp,'%3.0f'),'% accepted']);
    %     textloc(dtext,'northeast','FontSize',fontsz,'FontWeight',fontwg)
    %     xlabel('Time BP (a)');ylabel('\Delta T (K)');
    %     legend('Tikhonov/UPRE','MCMC (avg)','Luterbacher 2004','Heikkilae 2003','GSM ensemble','location','southeast')
    %     filename=strcat(tit,'_priorproxy');
    %     saveas(gcf,filename,plotfmt)ply, Reply
    %
    %--------------------------------------------------------------------------
end

