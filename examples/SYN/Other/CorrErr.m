clear all
close all
clc

pltpath='./';
datpath='./';
srcpath='../../';
utlpath='../../';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

N = 100;

L=7;
Err = 1.;
Cov=CovarGauss(Err*ones(N,1),L);

C    = chol(Cov);


err_nor =  Err*randn(N,1);
err_cor =  err_nor'*C;

figure 
plot(1:N,[err_nor(:)], '-b'); hold on
plot(1:N,[err_cor(:)], '-r','LineWidth',2); old on
