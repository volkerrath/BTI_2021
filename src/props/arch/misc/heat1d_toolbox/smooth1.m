function [s]=smooth1(f)
% SMOOTH1(f) smoothes input time series f with triangular  weights
% (1/4 1/2 1/4)
% V. R.  May. 28, 2003

[n1 n2]=size(f);

if n1==1,k=n2;else k=n1;end;
 

e=ones(k+4,1);
Op = spdiags([.25*e,.5*e,.25*e],-1:1,k+4,k+4);
p1=f(1);pk=f(k);

if n1==1, 
    fp=[p1 p1 f pk pk];
    sp=Op*fp';
else 
    fp=[p1;p1; f; pk; pk];
    sp=Op*fp;
end;
    
s=sp(3:k+2);

