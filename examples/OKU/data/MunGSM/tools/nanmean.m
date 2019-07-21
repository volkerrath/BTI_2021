function [avg,med]=nanmean(m)
% knanmean(X,[dim]) 
% returns the average and median of the columns for matrix X on ,
% treating NaNs as missing values. 
% VR july 2013
[i1 i2]=size(m);
a=isfinite(m);

avg=NaN(i1,1);

for col=1:i2
 v=m(a(:,col),col);
 avg(col)=mean(v);
 med(col)=median(v);
end
