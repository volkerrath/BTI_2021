clear all
close all
clc

% SET PATHS
addpath(['./tools'])

yr2s=31557600;s2yr=1/yr2s;

plots=0;
plotfmt='png';
linwdt=2;
fontwg='normal';
fontsz=16;
colscl='default';



disp('SET 2 - 2015  ')
nsamp=18;
sites={'SUL', 'VAT', 'ONE', 'OLK', 'FOR', 'VID', 'LAX'};
nsites=length(sites);whichsites=[1:nsites];
for isites=1:nsites
    site=sites{whichsites(isites)};
    nfilename{isites}=strcat([site,'_MunGSM_Data']);
end

TMP=dir('idBTsites.nn9*');
files={TMP.name};
whichfiles=[1:length(files)];nfiles=length(whichfiles);


for ifiles=1:nfiles
    file=files{whichfiles(ifiles)};
    disp(strcat([' Ensemble read from ' file ]));
    TMP=importdata(file, ' ', 1);
    TMP=TMP.data;
    t=TMP(:,1)*1000*yr2s;
    for isites=1:nsites
        site=sites{whichsites(isites)};
        offset=1+0*nsamp+(whichsites(isites));
        T=TMP(:,offset);
        offset=1+1*nsamp+(whichsites(isites));
        H=TMP(:,offset);
        offset=1+2*nsamp+(whichsites(isites));
        W=TMP(:,offset);
        nvarname=genvarname(site);
        if ifiles==1, eval([ nvarname '= [t];']);end
        eval([ nvarname '= [' nvarname ' T H W  ];'])
    end
end
disp('  ')
for isites=1:nsites
    site=sites{whichsites(isites)};
    nvarname=genvarname(site);
    eval(['DAT=' nvarname ';']);[n1,n2]=size(DAT);
    save(nfilename{isites},'DAT');
    disp(strcat([' Ensemble written to ' nfilename{isites} ]));
    tx=abs(DAT(:,1)/(1000*yr2s));
    Tx=DAT(:,2:3:n2);Tavg=mean(Tx,2);Tmed=median(Tx,2);
    Hx=DAT(:,3:3:n2);Havg=mean(Hx,2);Hmed=median(Hx,2);
    Wx=DAT(:,4:3:n2);Wavg=mean(Wx,2);Wmed=median(Wx,2);
    F= strcat([nvarname,'_Ensemble']);
    save (F,'tx','Tx','Tavg','Tmed','Hx','Havg','Hmed','Wx','Wavg','Wmed');
    
    
    if plots
        
        figure
        plot(tx,Tx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('\Delta T (K)');
        
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_T_Ensemble']);
        saveas(gcf,F,plotfmt)
        
        figure
        plot(tx,Hx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('Ice thickness (m)');
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_H_Ensemble']);
        saveas(gcf,F,plotfmt)
        figure
        plot(tx,Wx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('Water Index (-)');
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_W_Ensemble']);
        saveas(gcf,F,plotfmt)
        
    end
end
disp('  ')
disp('  ')
disp('SET 1 - 2013  ')

nsamp=9;
sites  ={'CO1', 'CO2', 'GRV', 'SG3','OKU', 'CX1', 'UDR', 'CZE', 'TOR'};
nsites=length(sites);whichsites=[1:nsites];
for isites=1:nsites
    site=sites{whichsites(isites)};
    nfilename{isites}=strcat([site,'_MunGSM_Data']);
end


TMP=dir('idBTsites.nn4*');
files={TMP.name};
whichfiles=[1:length(files)];nfiles=length(whichfiles);


for ifiles=1:nfiles
    file=files{whichfiles(ifiles)};
    disp(strcat([' Ensemble read from ' file ]));
    TMP=importdata(file, ' ', 1);
    TMP=TMP.data;
    t=TMP(:,1)*1000*yr2s;
    for isites=1:nsites
        site=sites{whichsites(isites)};
        offset=1+0*nsamp+(whichsites(isites));
        T=TMP(:,offset);
        offset=1+1*nsamp+(whichsites(isites));
        H=TMP(:,offset);
        offset=1+2*nsamp+(whichsites(isites));
        W=TMP(:,offset);
        nvarname=genvarname(site);
        if ifiles==1, eval([ nvarname '= [t];']);end
        eval([ nvarname '= [' nvarname ' T H W];'])
    end
end
disp('  ')

for isites=1:nsites
    site=sites{whichsites(isites)};
    nvarname=genvarname(site);
    eval(['DAT=' nvarname ';']);[n1,n2]=size(DAT);
    save(nfilename{isites},'DAT');
    disp(strcat([' Ensemble written to ' nfilename{isites} ]));
    tx=abs(DAT(:,1)/(1000*yr2s));
    Tx=DAT(:,2:3:n2);Tavg=mean(Tx,2);Tmed=median(Tx,2);
    Hx=DAT(:,3:3:n2);Havg=mean(Hx,2);Hmed=median(Hx,2);
    Wx=DAT(:,4:3:n2);Wavg=mean(Wx,2);Wmed=median(Wx,2);
    F= strcat([nvarname,'_Ensemble']);
    save (F, 'tx','Tx','Tavg','Tmed','Hx','Havg','Hmed','Wx','Wavg','Wmed');
    
    
    if plots
        figure
        plot(tx,Tx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('\Delta T (K)');
        
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_T_Ensemble']);
        saveas(gcf,F,plotfmt)
        
        figure
        plot(tx,Hx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('Ice thickness (m)');
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_H_Ensemble']);
        saveas(gcf,F,plotfmt)
        figure
        plot(tx,Wx,'LineWidth',linwdt);hold on
        grid on
        set(gca,'xscale','log','xdir','rev',...
            'FontSize',fontsz,'FontWeight',fontwg);
        xlabel('Time BP/2000 (ka)');
        ylabel('Water Index (-)');
        ptext= strcat(['Site : ',nvarname]);
        textloc(ptext,'north','FontSize',fontsz-1,'FontWeight',fontwg);
        F= strcat([nvarname,'_W_Ensemble']);
        saveas(gcf,F,plotfmt)
        
    end
end



%                 'xtick',[10 100 1000 10000 100000],...
%             'ytick',[-4 0 4 8],...
%               dtit=strrep(tit,'_',' ');

