clear all
close all
clc

T=[-8:0.1:2];

Tf=0.;w=1.;
[Theta1,dTheta1]=ftheta_general(T,Tf,w); 

Tf=0.;
[Theta2,dTheta2]=ftheta_galushkin(T,Tf); 

Ts=-0.2;b=.6;
[Theta3,dTheta3]=ftheta_nicolsky(T,Ts,b); 

ThetaAll=[Theta1' Theta2' Theta3'];

dThetaAll=[dTheta1' dTheta2' dTheta3'];


figure;
plot(T, ThetaAll,'LineWidth',2);grid on;
set(gca,'FontSize',14,'Fontweight', 'bold')
legend('Lunardini','Galushkin','Nicolsky','Location','NorthWest')
title('partition functions','FontSize',14);
xlabel('T (C)');
ylabel('fluid content (-)')
file='Test_ftheta.ps';saveas(gcf,file,'epsc2')


figure;
plot(T, dThetaAll,'LineWidth',2);grid on;
set(gca,'FontSize',14,'Fontweight', 'bold')
legend('Lunardini','Galushkin','Nicolsky','Location','NorthWest')
title('gradient partition functions','FontSize',14);
xlabel('T (C)');
ylabel('\delta fluid content (-)')
file='Test_dftheta.ps';saveas(gcf,file,'epsc2')