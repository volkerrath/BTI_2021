function [m]=wmean(v,w,mode)
% WMEAN weighted mean
% 
% [m]=wmean(v,w,mode) calculates the weighted arithmetic, 
% geometric, square-root, or harmonic mean of the rows of matrix v with
% weights w 
% For default  w = [1 ...] and mode='arithmetic'.
% V.R., A. H. May 15, 2003

if nargin<3, mode='arithmetic'; end
if nargin<2, [n1,n2]=size(v);w=ones(n1,n2); end

switch lower(mode)
    case 'arithmetic'
        m=sum(v.*w,2)./sum(w,2);
    case 'geometric'
        m=exp(sum(log(v).*w,2)./sum(w,2));
    case 'harmonic'
        m=1./(sum(w./v,2)./sum(w,2));
    case 'square-root'
        m = sum(w.*sqrt(v),2).^2./sum(w,2);
    otherwise
        disp(['WMEAN: mode set to arithmetic, >', mode, '<  not defined'])
        m=sum(v.*w,2)./sum(w,2);
        
end

