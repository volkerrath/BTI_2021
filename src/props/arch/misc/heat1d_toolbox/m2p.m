function [p]=m2p(m,flag)
% transforms from inverse parameter space m to physical parameter p
% V. R.,Nov 4, 2002
 p=m;
 pnt=find(flag==1);
 p(pnt) = exp(m(pnt)); 
