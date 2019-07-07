function [rci]=rciT(T,P)
% set pure water heat capacity
% VR RWTH Aachen Univerity,   April 25, 2004
T=T(:);
rci(T>0.001) = 1000*...
         (1.94+ ...
          7.14*T);
