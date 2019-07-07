function [m]=p2m(p,flag)
% transforms from physical parameter space p to inverse parameter m
% V. R.,Nov 4, 2002
 m=p;
 pnt=find(flag==1);
 m(pnt) = log(p(pnt)); 
