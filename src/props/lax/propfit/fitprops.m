clear all
close all
clc
%GRAPHICS
plotfmt='epsc2';
outsteps=false;
graph=false;
fontwg='normal';
fontsz=16;

Tfit=[0:1:120];


% OBSERVATIONS
Tobs = [25       60      100]; 
kobs = [2.91    2.82    2.72];
robs = [2743    2743   2743];
cobs = [712      764    824];

% THERMAL CONDUCTIVITY
K=[0.7 770];
Kfit = nlinfit(Tobs,kobs,@kfun,K);
kfit=kfun(Kfit,Tfit);

disp([' ']);
disp(['Thermal conductivity (Clauser & Huenges, 1995): ']);
disp(['Original coeffcient values = ',num2str(K   ,' %7.3g  %7.3g')]);
disp(['Fitted coeffcient values   = ',num2str(Kfit,' %7.3g  %7.3g')]);


figure
plot(Tobs,kobs,'+b','LineWidth',2); hold on
plot(Tfit,kfit,'-r','LineWidth',2); hold on
set(gca,'FontSize',fontsz,'FontWeight',fontwg);
grid on
tickstr={'measured', 'fit'};legend(tickstr,'location','southwest');
ylabel('\lambda (W m^{-1}K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
filename=strcat(['Olkiluoto_TC']);
saveas(gcf,filename,plotfmt)

% HEAT CAPACITY
C = [0.0044 2.134];
Cfit = nlinfit(Tobs,cobs,@cfun,C);
cfit=cfun(Cfit,Tfit);
disp([' ']);
disp(['Heat capacity (Mottaghgy et al. 2005): ']);
disp(['Original coeffcient values = ',num2str(C   ,' %7.3g  %7.3g')]);
disp(['Fitted coeffcient values   = ',num2str(Cfit,' %7.3g  %7.3g')]);


figure
plot(Tobs,cobs,'+b','LineWidth',2); hold on
plot(Tfit,cfit,'-r','LineWidth',2); hold on
grid on;
set(gca,'FontSize',fontsz,'FontWeight',fontwg);
tickstr={'measured', 'fit'};legend(tickstr,'location','southeast');
ylabel('c_p (J kg^{-1} K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg);
filename=strcat(['Olkiluoto_CP']);
saveas(gcf,filename,plotfmt)



