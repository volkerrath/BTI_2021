function [vc] = n2c(vn,d)
% VC=N2C(VN) interpolates node-to-cell centers
%
% function [vc] = n2c(vn,z) interpolates nodewise
% defined parameter to cells. it takes node values vn and
% cell sizes d, gives cell-centered values vc. size of vc is
% length(vn)-1=length(d).
% V. R.,  Oct. 24, 2001 

nn=length(vn);nc=nn-1;    
vc =0.5*(vn(1:nc)+vn(2:nc+1));
