function [cpmT]=cpmT(cp,T,T0,T1)
% calculates heat capacity of rocks as function of temperature
% based on the formula given in Cermak & Rybach (Landolt-Boernstein, 
% Vol. V/1a Physical Properties of rocks) 
% T           =  temperature in C  acording to flag unit
% vr april 25, 2004
A=754.;B=6.14e-4;C=1.928e4;
if nargin < 3, T0=20.;T1=800; end
if nargin < 3, T0=20.;T1=800; end

[n1 n2]=size(T);if n1==1, T=T'; end 
[n1 n2]=size(cp);if n1==1, cp=cp'; end

 
TK=T+273.15;T0K=T0+273.15;T1K=T1+273.15;
% cpmT=ones(size(cp))*800;
cp0=A*(1+B*T0K-C./T0K^2);f=cp./cp0;
cp1=A*(1+B*T1K-C./T1K.^2);

cp=A*(1+B*TK-C./(TK.*TK));
ilimit=find(T>=T1);cp(ilimit)=cp1;
% cpmT=cp;
cpmT=f.*cp;

