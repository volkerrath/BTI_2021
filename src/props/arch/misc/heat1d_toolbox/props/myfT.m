function [my,T]=myfT1(T)
% VR RWTH Aachen University,   April 25, 2004

[n1 n2]=size(T);if n1==1, T=T'; end 
% Sutra:    A =239.4e-7;    B=248.37;   C=133.16;
% HST:      A =243.18e-7;   B=247.8;    C=133.16;
%T=T+133.16; % T in C 
%T=T-140.;   % T in K 
A =243.18e-7;B=247.8;C=133.16;
my=A*10.^(B./(T+C));
