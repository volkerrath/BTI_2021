% test we1d
clear aall; close all;
opts = struct('Format','psc2','Color','RGB','Resolution',600);


M=100;
K=M+1;

dz=1.;
zn=dz*[0:1:M];
m = ones(M,1); 
m(21:40)=3; 
m(41:65)=1;
m(66:100)=2;

% figure
% stairs(z(1:M),m); grid on;
% set(gca,'XScale','log')
% xlabel('x');ylabel('m');ylim([-1 4]);
% title(['true model'])
% filename='Truemodel.ps'
% exportfig(gcf,filename,opts)





eps=1.e-6;
bet=1.e-6;
zm=0.5*(zn(2:K)+zn(1:M));

[W1]=we1d(m,zm,'SM',eps,bet) ;ww1=spdiags(W1);
[W2]=we1d(m,zm,'ME',eps,bet) ;ww2=spdiags(W2);
[W3]=we1d(m,zm,'TV',eps,bet) ;ww3=spdiags(W3);
[W4]=we1d(m,zm,'MS',eps,bet) ;ww4=spdiags(W4);
[W5]=we1d(m,zm,'MGS',eps,bet);ww5=spdiags(W5);

% [W0,W0i]=we1d(m,zm,'NIX',eps,bet);ww0=diag(W0);
% [W1,W1i]=we1d(m,zm,'SM',eps,bet) ;ww1=diag(W1);
% [W2,W2i]=we1d(m,zm,'ME',eps,bet) ;ww2=diag(W2);
% [W3,W3i]=we1d(m,zm,'TV',eps,bet) ;ww3=diag(W3);
% [W4,W4i]=we1d(m,zm,'MS',eps,bet) ;ww4=diag(W4);
% [W5,W5i]=we1d(m,zm,'MGS',eps,bet);ww5=diag(W5);

figure
plot(zm,[ww1 ww2 ww3 ww4 ww5]); 
legend('SM','ME','TV', 'MS', 'MGS')
% set(gca,'XScale','log')
xlabel('x');ylabel('W_e');
% ylim([-1 4]);
title(['weights']);
filename='Weights.ps';
exportfig(gcf,filename,opts)
