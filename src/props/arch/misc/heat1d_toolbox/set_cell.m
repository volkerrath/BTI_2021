function [l]=set_cell(zl,lam,z,method) 
if nargin < 4, method='nearest'; end
lz=length(zl);lg=length(z);

zm(1:lg-1)=0.5*(z(1:lg-1)+z(2:lg));lm(1:lg-1)=NaN;

for i=1:lg-1
    in=find(zl>z(i) & zl <= z(i+1));
    if ~isempty(in)
        n=length(in);lm(i)=sum(lam(in))/n;
    end       
end 

points=(find(~isnan(lm)));

l=interp1(zm(points),lm(points),zm,method,'extrap');



