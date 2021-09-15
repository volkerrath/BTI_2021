function [rcf]=rcfT(T,P)
% set pure water heat capacity
% VR RWTH Aachen Univerity,   April 25, 2004
T=T(:);
rcf=rhofT(T).*cpfT(T);