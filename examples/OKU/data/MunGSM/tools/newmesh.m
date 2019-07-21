clear all
close all
clc

load OLK_obs

z1=[0 zT'];dz1=diff(z1);
dz2=2.5*1.03.^[1:128];

dz=[dz1 dz2]; z=cumsum(dz)

save('OLK_mesh','z','dz')


 