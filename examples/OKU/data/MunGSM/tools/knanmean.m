function res=knanmean(m,dim)
% knanmean(X,[dim]) 
% returns the average of n-dimension matrix X on dimension dim (default=1),
% treating NaNs as missing values. 
% Alex 2009

if nargin==1, dim=1; end

a=isnan(m);
co=sum(~a,dim);

m(a)=0;
s=sum(m,dim);
res=s./co;
res(isinf(res))=nan;