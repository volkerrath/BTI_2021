function [cpmT]=cpmT(cpm0,T)
% calculates heat capacity of rocks as function of temperature
% based on the formula given in Cermak & Rybach (Landolt-Boernstein,
% Vol. V/1a Physical Properties of rocks)
% T           =  temperature in C  acording to flag unit
% vr april 25, 2004
A=754.;B=6.14e-4;C=1.928e4;
T0=20.;T1=800;

T=T(:);cpm0=cpm0(:);

TK=T+273.15;T0K=T0+273.15;T1K=T1+273.15;
% cpmT=ones(size(cp))*800;
cp0=A*(1+B*T0K-C./T0K^2);f=cpm0./cp0;
cp1=A*(1+B*T1K-C./T1K.^2);

cp=A*(1+B*TK-C./(TK.*TK));
ilimit=find(T>=T1);cp(ilimit)=cp1;
% cpmT=cp;
cpmT=f.*cp;
cpmT=cpmT(:);

