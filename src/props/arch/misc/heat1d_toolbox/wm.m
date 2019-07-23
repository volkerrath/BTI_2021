function [W]=wm(J)
% [W]=function wm(J) calculates the diagonal weighting 
% matrix for the parameters
% Wm
% =>     J          Jacobian (Nd x Np)
% <=     W          weight sparse diagonal matrix
% see:
% Zhdanov, M. S. (2002): Geophysical inverse theory and 
%       regularization problems, p. 80ff 

[Nd Np]=size(J);
A=sum(J.^2);W=sqrt(sqrt(A))';
W=spdiags(W,[0],Np,Np);
