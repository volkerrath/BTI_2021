clear all, close all;

T=[1:1:150];

figure;
[c]=cpt(1,T,20);
plot(T,c);hold on;

title('c_{p} = f(T) (Cermak & Rybach, 1982)');grid on;
xlabel('T (C)');ylabel('c_{p} (J K^{-1} kg^{-1})')
opts = struct('Format','psc2','Color','rgb');
exportfig(gcf,'cp_cermak.eps',opts)
