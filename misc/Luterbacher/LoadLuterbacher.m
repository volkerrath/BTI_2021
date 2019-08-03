clear all
close all
clc
addpath([pwd,filesep,'tools']);
exportToZip(mfilename,strcat([mfilename '_' date '.zip']));

plot_fmt='png';
fontwg='normal';
fontsz=14;
plot_4=0;
plot_s=0;


nlon=130;
longitude =[...
   -24.75, -24.25, -23.75, -23.25, -22.75,  -22.25, -21.75, -21.25, -20.75, -20.25, ...
   -19.75, -19.25, -18.75, -18.25, -17.75,  -17.25, -16.75, -16.25, -15.75, -15.25, ...
   -14.75, -14.25, -13.75, -13.25, -12.75,  -12.25, -11.75, -11.25, -10.75, -10.25, ...
    -9.75,  -9.25 , -8.75,  -8.25,  -7.75,  -7.25,  -6.75,  -6.25,  -5.75,  -5.25, ...
    -4.75,  -4.25,  -3.75,  -3.25,  -2.75,  -2.25,  -1.75,  -1.25,  -0.75,  -0.25, ...
     0.25,   0.75,   1.25,   1.75,   2.25,   2.75,   3.25,   3.75,   4.25,   4.75, ...
     5.25,   5.75,   6.25,   6.75,   7.25,   7.75,   8.25,   8.75,   9.25,   9.75, ...
    10.25,  10.75,  11.25,  11.75,  12.25,  12.75,  13.25,  13.75,  14.25,  14.75, ...
    15.25,  15.75,  16.25,  16.75,  17.25,  17.75,  18.25,  18.75,  19.25,  19.75, ...
    20.25,  20.75,  21.25,  21.75,  22.25,  22.75,  23.25,  23.75,  24.25,  24.75, ...
    25.25,  25.75,  26.25,  26.75,  27.25,  27.75,  28.25,  28.75,  29.25,  29.75, ...
    30.25,  30.75,  31.25,  31.75,  32.25,  32.75,  33.25,  33.75,  34.25,  34.75, ...
    35.25,  35.75,  36.25,  36.75,  37.25,  37.75,  38.25,  38.75,  39.25,  39.75  ...
   ];
nlat=70;
latitude = [...
   35.25, 35.75, 36.25, 36.75, 37.25, 37.75, 38.25, 38.75, 39.25, 39.75, ...
   40.25, 40.75, 41.25, 41.75, 42.25, 42.75, 43.25, 43.75, 44.25, 44.75, ...
   45.25, 45.75, 46.25, 46.75, 47.25, 47.75, 48.25, 48.75, 49.25, 49.75, ...
   50.25, 50.75, 51.25, 51.75, 52.25, 52.75, 53.25, 53.75, 54.25, 54.75, ...
   55.25, 55.75, 56.25, 56.75, 57.25, 57.75, 58.25, 58.75, 59.25, 59.75, ...
   60.25, 60.75, 61.25, 61.75, 62.25, 62.75, 63.25, 63.75, 64.25, 64.75, ...
   65.25, 65.75, 66.25, 66.75, 67.25, 67.75, 68.25, 68.75, 69.25, 69.75  ...
   ];



load list.dat
[nl,nv]=size(list);nval=list(1,4);

load temp

[nrec,ndum]=size(temp);
Tseas=reshape(temp,nlat,nrec/nlat,nlon);
Tseas(Tseas>=100)=NaN;
[nlat2,nseas,nlon2]=size(Tseas);
Tseas=permute(Tseas,[1 3 2]);
size(Tseas);
lseas(:,1)=floor(list(:,1)/10000);
lseas(:,2)=list(:,1)-lseas*10000;
filename='Luterbacher2004_seas.mat';
save(filename,'Tseas','nseas','lseas')

ycount=1;
for yr=1:4:nseas-1
    Tyear(:,:,ycount) =knanmean(Tseas(:,:,yr:yr+3),3);
    ycount=ycount+1;
end

[nyear,nlon3,nlat3]=size(Tyear);
lyear=floor(list(1:4:nseas,1)/10000);
filename='Luterbacher2004_year.mat';
save(filename,'Tyear','nyear','lyear')
lyear=lyear+0.5;

% OKU coordinates
SITE='OKU';
OlatG=62; OlatM=43; OlatS=4;
OlonG=29; OlonM=3;  OlonS=43;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_OKU = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_OKU_10 = slidefun (@mean, W, T_OKU, 'central'); 
W=30;
T_OKU_30 = slidefun (@mean, W, T_OKU, 'central');




YL=lyear;TL=T_OKU;T10L=T_OKU_10;T30L=T_OKU_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_OKU,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_OKU_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_OKU_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['OKU  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['OKU_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_OKU_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['OKU  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['OKU_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end


% KOLA 
SITE='SG3';
OlatG=69; OlatM=23; OlatS=0;
OlonG=30; OlonM=36; OlonS=0;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_SG3 = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_SG3_10 = slidefun (@mean, W, T_SG3, 'central'); 
W=30;
T_SG3_30 = slidefun (@mean, W, T_SG3, 'central');



YL=lyear;TL=T_SG3;T10L=T_SG3_10;T30L=T_SG3_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_SG3,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_SG3_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_SG3_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['SG3  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['SG3_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_SG3_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['SG3  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['SG3_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end


% UDR coordinates: % 54?14?N 22?56?E
SITE='UDR';
OlatG=54; OlatM=14; OlatS=0;
OlonG=22; OlonM=56; OlonS=0;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_UDR = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_UDR_10 = slidefun (@mean, W, T_UDR, 'central'); 
W=30;
T_UDR_30 = slidefun (@mean, W, T_UDR, 'central');

YL=lyear;TL=T_UDR;T10L=T_UDR_10;T30L=T_UDR_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_UDR,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_UDR_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_UDR_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['UDR  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['UDR_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_UDR_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    %ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['UDR  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['UDR_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end





% CZE coordinates:  52 6 0N     17 30 0E
SITE='CZE';
OlatG=52; OlatM=06; OlatS=0;
OlonG=17; OlonM=30; OlonS=0;
% convert to decimal grades
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2);
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);
T21=Tyear(lap,lom,:);
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_CZE = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_CZE_10 = slidefun (@mean, W, T_CZE, 'central');
W=30;
T_CZE_30 = slidefun (@mean, W, T_CZE, 'central');

YL=lyear;TL=T_CZE;T10L=T_CZE_10;T30L=T_CZE_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_CZE,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CZE_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_CZE_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CZE  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CZE_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CZE_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    %ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CZE  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CZE_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end
 
% TOR coordinates:     53 0 0N  18 48 0E
SITE='TOR';
OlatG=54; OlatM=14; OlatS=0;
OlonG=18; OlonM=48; OlonS=0;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_TOR = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_TOR_10 = slidefun (@mean, W, T_TOR, 'central'); 
W=30;
T_TOR_30 = slidefun (@mean, W, T_TOR, 'central');

YL=lyear;TL=T_TOR;T10L=T_TOR_10;T30L=T_TOR_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_TOR,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_TOR_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_TOR_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['TOR  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['TOR_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_TOR_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    %ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['TOR  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['TOR_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end
 

% CO1 coordinates:      63 23 60N 13 11 60"E
SITE='CO1';
OlatG=63; OlatM=24; OlatS=0;
OlonG=13; OlonM=12; OlonS=0;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_CO1 = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_CO1_10 = slidefun (@mean, W, T_CO1, 'central'); 
W=30;
T_CO1_30 = slidefun (@mean, W, T_CO1, 'central');

YL=lyear;TL=T_CO1;T10L=T_CO1_10;T30L=T_CO1_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_CO1,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO1_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO1_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CO1  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CO1_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO1_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    %ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CO1  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CO1_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end
 

% CO2 coordinates:     63°18 0 N   13 30 0E
SITE='CO2';
OlatG=63; OlatM=18; OlatS=0;
OlonG=13; OlonM=30; OlonS=0;
% convert to decimal grades                                             
mla = OlatM / 60; sla = ( OlatS / 60 ) / 60;la = OlatG + mla + sla;
mlo = OlonM / 60; slo = ( OlonS / 60 ) / 60;lo = OlonG + mlo + slo;

[x,y,utmzone] = deg2utm(la,lo);   

% determine cell for interpolation
lop=min(find(longitude>lo));lom=lop-1;Glop=longitude(lop);Glom=longitude(lom);
lap=min(find(latitude >la));lam=lap-1;Glap=latitude(lap);Glam=latitude(lam);

[x11,y11] = deg2utm(lam,lom); r11=1./(x11^2+y11^2); 
[x12,y12] = deg2utm(lap,lom); r12=1./(x12^2+y12^2);  
[x21,y21] = deg2utm(lam,lop); r21=1./(x21^2+y21^2);  
[x22,y22] = deg2utm(lap,lop); r22=1./(x22^2+y22^2);
scl=r11+r12+r21+r22;
% extract cell values 
T11=Tyear(lam,lom,:);
T22=Tyear(lap,lop,:);
T12=Tyear(lam,lop,:);  
T21=Tyear(lap,lom,:);  
T11=T11(:);T12=T12(:);T21=T21(:);T22=T22(:);
T_CO2 = (T11*r11+T12*r12+T21*r21+T22*r22)/scl;

W=10;
T_CO2_10 = slidefun (@mean, W, T_CO2, 'central'); 
W=30;
T_CO2_30 = slidefun (@mean, W, T_CO2, 'central');

YL=lyear;TL=T_CO2;T10L=T_CO2_10;T30L=T_CO2_30;
filename=strcat(['Luterbacher2004_',SITE]);
save(filename,'YL','TL','T10L','T30L')
message=strcat(['Recons for ',SITE, ' written to ',filename]);disp(message)

if plot_s
    figure
    plot(lyear,T_CO2,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO2_10,'Color','b','LineWidth',2,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO2_30,'Color','r','LineWidth',2,'LineStyle','-');hold on;grid on
    legend('orig','avg 10 yrs','avg 30 yrs','Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CO2  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CO2_LuterbacherYear']);
    saveas(gcf,file,plot_fmt)
end

if plot_4
    figure
    plot(lyear,T11,'Color','k','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T21,'Color','r','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T12,'Color','g','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T22,'Color','b','LineWidth',1,'LineStyle','-');hold on;grid on
    plot(lyear,T_CO2_30,'Color','m','LineWidth',3,'LineStyle','-');hold on;grid on
    l1=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glom)]);
    l2=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glom)]);
    l3=strcat(['lat= ',num2str(Glam),' lon= ',num2str(Glop)]);
    l4=strcat(['lat= ',num2str(Glap),' lon= ',num2str(Glop)]);
    legend(l1,l2,l3,l4,'Location','southwest')
    
    xlabel('year A.D.','FontSize',fontsz,'FontWeight',fontwg);
    ylabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
    xlim([1490 2010]);
    %ylim([-4 6])
    
    set(gca,'FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['CO2  (from Luterbacher, 2004)']),'FontSize',fontsz,'FontWeight',fontwg)
    
    file=strcat(['CO2_LuterbacherYear_Cell']);
    saveas(gcf,file,plot_fmt)
end
 


