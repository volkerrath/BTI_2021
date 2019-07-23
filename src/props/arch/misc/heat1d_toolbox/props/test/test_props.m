clear all;close all;
opts = struct('Format','psc2','Color','rgb','Resolution',600);

load('IAP.dat')


z=[0:10:3000];

P0=0.01;
P=P0+998.*9.81*z;

T0=0.01;
T=T0+0.03*z;

%figure;
%subplot(1,2,1)
%plot(T,z,'-r');grid on;
%xlabel('T (C)');ylabel('z (m)');
%set(gca,'ydir','reverse')
%title('temperature ');
%subplot(1,2,2)
%plot(P*1.e-6,z,'-r');grid on;
%xlabel('P (MPa)');ylabel('z (m)');
%set(gca,'ydir','reverse')
%title('pressure');
%filename='Props_1.ps';
%exportfig(gcf,filename,opts);%close(gcf)


% rhof

r1=rhofT_batzle(T,P)';
r2=rhofT_fehmn(T,P)';
r3=rhofT_mavko(T,P)';
r4=IAP(:,4);
figure;
plot([r1],z,'-r','linewidth',1);grid on;hold on;
plot([r2],z,'-g','linewidth',1);grid on;hold on;
plot([r3],z,'-b','linewidth',1);grid on;hold on;
plot([r4],z,'-m','linewidth',2);grid on;hold on;
set(gca,'ydir','reverse')
title('density'); % xlim([960 1010])
xlabel('(\rho (kg/m^3)');ylabel('z (m)');
legend('Batzle 1992','FEHMN 1994','Mavko 1998','IAP 1995',2)
filename='Rhof.ps';exportfig(gcf,filename,opts);%close(gcf)

% lambdaf
r1=kfT_phillips(T)';
r2=kfT_ramirez(T)';
r3=IAP(:,5)
figure;
plot([r1,],z,'-r','linewidth',1);grid on;hold on;
plot([r2,],z,'-b','linewidth',1);grid on;hold on;
plot([r3,],z,'-m','linewidth',2);grid on;hold on;
%plot([r3],z);grid on;
set(gca,'ydir','reverse')
title('thermal conductivity'); % xlim([960 1010])
xlabel('(\lambda (W/mK)');ylabel('z (m)');
legend('Phillips 1981','Ramirez 1994','IAP 1995',1)
filename='Lamf.ps';exportfig(gcf,filename,opts);%close(gcf)

% cpf
r1=cpfT_fehmn(T,P)';
r3=IAP(:,6)

figure;
plot([r1,],z,'-b','linewidth',1);grid on;hold on;
plot([r3,],z,'-m','linewidth',2);grid on;hold on;
set(gca,'ydir','reverse')
title('isobaric heat capacity'); % xlim([960 1010])
xlabel('(c_p (J/kg)');ylabel('z (m)');
legend('FEHMN 1994','IAP 1995',2)
filename='Cpf.ps';exportfig(gcf,filename,opts);%close(gcf)
