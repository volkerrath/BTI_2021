function [S]=threshsp(F,thresh)
% THRESHSP contracts full (covariance) matrix to sparse matrix S
% given threshold thresh
% V. R. Oct. 15, 2001 

if nargin < 2, thresh=0.001; end

is=find(abs(F)>thresh);v=F(is);
[i,j]=ind2sub(size(F),is);
S=sparse(i,j,v);
