%GRAPHICS
plotfmt='epsc2';
linwdt=2;
fontwg='normal';
fontsz=16;

addpath('./src')


phi=[0:0.01:0.5];
kf=.65*ones(size(phi));
ke=1.8*ones(size(phi));

kma=ke2km(ke,kf,phi,'a')
kms=ke2km(ke,kf,phi,'s')
kmg=ke2km(ke,kf,phi,'g')
kmh=ke2km(ke,kf,phi,'h')



figure
plot(phi,[kma(:),kms(:),kmg(:),kmg(:)],'LineWidth',linwdt)
grid on
set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
legend(' arithmetic','square root','geometric','harmonic','location','northwest')
xlabel('porosity (-)','FontSize',fontsz,'FontWeight',fontwg);
ylabel('\lambda (W m^{-1}K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
title('\lambda_m from mixing laws');
file=strcat(['LambdaMix']);
saveas(gcf,file,plotfmt)
