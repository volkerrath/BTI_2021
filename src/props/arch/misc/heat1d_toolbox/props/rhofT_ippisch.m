function [rhof]=rhofT(T)
% calculates density of water as function of temperature
% based on the formula given in Ippisch(2001)
% T           =  temperature in C 
% vr april 19, 2004
[n1,n2]=size(T); if n1~=1, T=T'; end 

       B=3.17e-4;G=2.56e-6;rref=1000.3;Tref=0; 
       delT= T - Tref;
       rhof= 1000.3 * (1.0 - delT*B - delT.*delT*G);
%       rhof=1000.3*ones(size(T));
