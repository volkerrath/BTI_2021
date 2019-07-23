clear all
close all
clc

% GRAPHICS
fontwg='normal';
fontsz=16;
linwdt=2;



T=[-50:1:0];


k1=kiT_ling(T);
k2=kiT_fukosako(T);

figure;
plot(T,[k1 k2],'LineWidth',linwdt)
hold on; grid on;
titstrng='Ice Thermal Conductivity';title(titstrng,'FontSize',fontsz,'FontWeight',fontwg);
set(gca,'FontSize',fontsz,'FontWeight',fontwg)
ylabel('\lambda');xlabel('T (C)');
legend('std','Fukosako (1990) ',...
'location','southwest','orientation','horizontal')
file=strcat(['IceConductivity.eps']);
saveas(gcf,file,'epsc2')


r1=rhoiT(T);
r2=rhoiT_fukosako(T);
figure;
plot(T,[r1 r2],'LineWidth',linwdt)
hold on; grid on;
titstrng='Ice Density';title(titstrng,'FontSize',fontsz,'FontWeight',fontwg);
set(gca,'FontSize',fontsz,'FontWeight',fontwg)
ylabel('\rho (kg/m^3)');xlabel('T (C)');
legend('std','Fukosako (1990) ',...
'location','southwest','orientation','horizontal')
file=strcat(['IceDensity.eps']);
saveas(gcf,file,'epsc2')



c1=cpiT(T);
c2=cpiT_fukosako(T);
figure;
plot(T,[c1 c2],'LineWidth',linwdt)
hold on; grid on;
titstrng='Ice Heat Capacity';title(titstrng,'FontSize',fontsz,'FontWeight',fontwg);
set(gca,'FontSize',fontsz,'FontWeight',fontwg)
ylabel('c_p ()');xlabel('T (C)');
legend('std','Fukosako (1990) ',...
'location','southwest','orientation','horizontal')
file=strcat(['IceCp.eps']);
saveas(gcf,file,'epsc2')



