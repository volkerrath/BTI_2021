clear all;close all;clc;
opts = struct('Format','psc2','Color','rgb','Resolution',600);


z=[0:10:3000];

P0=0.01;
P=P0+998.*9.81*z;

T0=0.01;
T=T0+0.03*z;

figure;
subplot(1,2,1)
plot(T,z,'-r');grid on;
xlabel('T (C)');ylabel('z (m)');
set(gca,'ydir','reverse')
title('temperature ');
subplot(1,2,2)
plot(P*1.e-6,z,'-r');grid on;
xlabel('P (MPa)');ylabel('z (m)');
set(gca,'ydir','reverse')
title('pressure');
filename='Props_1.ps';
exportfig(gcf,filename,opts);%close(gcf)



% lambdaf
k0=3.;
r1=kmT_haenel(k0,T)';
r2=kmT_lehmann(k0,T)';
r3=kmT_sass(k0,T)';
r4=kmT_vost(k0,T)';

figure;
plot([r1,],z,'-r','linewidth',2);grid on;hold on;
plot([r2,],z,'-b','linewidth',2);grid on;hold on;
plot([r3,],z,'-m','linewidth',2);grid on;hold on;
plot([r4,],z,'-k','linewidth',2);grid on;hold on;
%plot([r3],z);grid on;
set(gca,'ydir','reverse')
title('thermal conductivity'); % xlim([960 1010])
xlabel('(\lambda (W/mK)');ylabel('z (m)');
legend('Haenel','Lehmann','Sass','vosteen','position','southeast')
filename='Lam.ps';exportfig(gcf,filename,opts);%close(gcf)

% cpm
cp0=750.
r1=cpmT_cermak(cp0,T)';
r3=cpmT_hermann(cp0,T)';
% 
figure;
plot([r1,],z,'-b','linewidth',1);grid on;hold on;
plot([r3,],z,'-m','linewidth',2);grid on;hold on;
set(gca,'ydir','reverse')
title('isobaric heat capacity'); % xlim([960 1010])
xlabel('(c_p (J/kg)');ylabel('z (m)');
legend('Cermak 1981','Hermann 1999',2)
filename='Cpm.ps';exportfig(gcf,filename,opts);%close(gcf)
