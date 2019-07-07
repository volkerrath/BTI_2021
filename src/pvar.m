function [p] = pvar(p0,pfac)
% PVAR calculates new test regularization parameters
if nargin <2 ,pfac=ones(size(p0)), end
p=p0(:).*pfac(:);
end

