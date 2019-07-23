function [cpf]=cpfT(T,P)
% set pure water heat capacity 
% VR RWTH Aachen Univerity,   April 25, 2004
[n1,n2]=size(P); if n1~=1, T=T'; end 
cpf=4217.6*ones(size(T));
