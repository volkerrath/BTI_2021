clear all 
close all
clc

syms  p Tb2 Tb1 T

f=p./(Tb1-Tb2);
Theta=(exp(-f*T)-exp(-f*Tb2))./(exp(-f*Tb1)-exp(-f*Tb2));

dTheta=diff(Theta,'T') 
simple(dTheta)
