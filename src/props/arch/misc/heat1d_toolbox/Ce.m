function [C]= Ce(n,s,L)
% calculates markovian (exponential) covariance C
s2=s*s;
C=zeros(n,n);
for i=1:n
    for j=i:n
        C(i,j)=exp(-abs(i-j)/L);
        C(j,i)=C(i,j);
    end
end
C=C*s2;
