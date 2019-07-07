function [cpf]=cpfT(T,P)
% set pure water heat capacity
% VR RWTH Aachen Univerity,   April 25, 2004
T=T(:);
cpf=4217.6*ones(size(T));
