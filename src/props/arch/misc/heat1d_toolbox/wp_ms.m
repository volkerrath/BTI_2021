function [Wp]=Wp_MS(m,ma,beta,method)
% calculates the minimum support (MS) weighting function
% for nonlinear IRLS 
%
% (Portiaguine & Zhdanov, 1999; Zhdanov, 2002)
% VR  May 15, 2004
l=length(m);
if nargin < 3, eps =1.e-5; method=1; end

e2=eps*eps;
switch method
   case {1}
       dm=m-ma;w=sqrt(dm.*dm+e2)';
   case {2}
       [L1] = l_1d(l,'l1');
       gm=L1*(m-ma)';w=sqrt(gm.*gm+e2);
end
Wp=spdiags(1./w,[0],l,l);


