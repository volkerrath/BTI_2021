clear all;
close all;

Np=50;
Nd=250;

J=rand(Nd,Np);
J(30:177,20:33)=J(30:177,20:33)*1000;
figure
imagesc(J)
colorbar

W=wm(J);
W1=inv(W);WW=W1*W1;

J2=J*WW;
figure
imagesc(J2)
colorbar
