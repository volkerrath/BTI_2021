clear all;close all;
p=[1:1:200]; par(1)=2; par(2)=3;
tic;[Cpp,sqCpp]=cov1_pp(p,par,'gauss');toc;
figure;imagesc(Cpp);colorbar;
tic;Cpp=threshsp(Cpp,.1);toc;
figure;imagesc(Cpp);colorbar;
tic;Cppi=inv(Cpp);toc;
figure;imagesc(Cppi);colorbar;
