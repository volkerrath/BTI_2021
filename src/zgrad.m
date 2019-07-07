    function [Qz,zq,dTdz]=Qz(T,z,k,smooth)
    % calculates  vertical gradient of T
    if nargin<4, smooth=0; end
    T=T(:); z=z(:);nz=length(z);
    zq=0.5*(z(1:nz-1)+z(2:nz));
    dTdz=diff(T(:))./diff(z(:));
    Qz  = -k.*dTdz;
    if smooth>0,
        spoints=linspace(min(zq),max(zq),smooth);
        pp=Qz(:);ss=splinefit(zq,pp,spoints,'r',0.5);
        Qz=ppval(ss,zq);
    end
%     Qz  = -k.*dTdz;
    end 
