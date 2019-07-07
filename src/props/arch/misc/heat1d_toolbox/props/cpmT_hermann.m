function [cpmT]=cpmT(cp,T,T0,T1,A,B,C)
% calculates heat capacity of rocks as function of temperature
% based on the formula given in Hermann(1999)
% T           =  temperature in C  acording to flag unit
% A           =  regression coeff., 0.7         < A < 0.8
% B           =  regression coeff., 0.0014      < B < 0.0022
% C           =  regression coeff., -0.00000334 < C < 0.0000016
% vr april 25, 2004
[n1 n2]=size(T);if n1==1, T=T'; end 
if nargin < 5, A=0.75e3;B=0.0018e3;C=-0.00000245e3; end
if nargin < 3, T0=20.;T1=800; end
[n1 n2]=size(cp);if n1==1, cp=cp'; end 
[n1 n2]=size(T);if n1==1, T=T'; end 

cp0=A+B*T0+C*T0^2;f=cp./cp0;
cpmT =A+B*T +C*T.*T;
cp1=A+B*T1+C*T1^2;ilimit=find(T>=T1);cpmT(ilimit)=cp1;

