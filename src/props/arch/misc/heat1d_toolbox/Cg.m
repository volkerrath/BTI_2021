function [C]= Ce(n,s,L)
% calculates gaussian covariance C
s2=s*s;
C=zeros(n,n);
L2=2*L^2;
for i=1:n
    for j=i:n
        C(i,j)=exp(-abs(i-j)^2/L2);
        C(j,i)=C(i,j);
    end
end
C=C*s2;