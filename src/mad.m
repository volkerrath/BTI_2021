function [MAD]=mad(val)
% median of absolute deviations.
% vr 02/03/2013 05:48:10 PM 
MAD=1.4826*median(abs(val-median(val)));
