function [td,zd] = heat1dan(k,h,Qb,Ts,index,dz,extra)
% HEAT1DAN solves analytically the stationary heat equation
% (general layered case)
%
% [td,qd] = heat1dan(par,extra) calculates analytically
% temperatures, given a layered model for thermal conductivity
% and heat production.
% par (1:np) is layer thermal conductivity, par(np+1:2*np) is heat
% production.par(2*np+1:2*np+1) is basal heat flows for the sites and
% par(2*np+nsites:2*np+2) is surface temperatures.
%
% V. R., Sept.4, 2001

nlp=length(par)-2;np=nlp/2;
nl=length(index);
lambda =k(index(1:nl));     % layer thermal conductivity
hprod  =h(index(1:nl));     % layer heat production

if nargin > 3
    nd=length(extra);
    zd=extra(1:nd);           % depths for temperature evaluation
end
% calculate heat flow for each layer boundary
% by upward integration from base value
ql(nl+1)=Qb;
for i=nl+1:-1:2
    ql(i-1)=ql(i)+h(i-1)*dz(i-1);
end

% now calculate  temperatures and depth of layer boundaries
tl(1)=Ts;zl(1)=0;
for i=2:nl+1
    tl(i)=tl(i-1)+ql(i-1).*dz(i-1)/lambda(i-1)-0.5*hprod(i-1)*dz(i-1)*dz(i-1)/lambda(i-1);
    zl(i)=zl(i-1)+dz(i-1);
end

% interpolate  temperatures  and heat flows

if nargin > 3
    for i=1:nd
        for j=1:nl
            if zd(i) > zl(j) &    zd(i) <= zl(j+1)
                dzd=zd(i)-zl(j);
                td(i)=tl(j)+ql(j)*dzd/lambda(j)-0.5*hprod(j)*dzd*dzd/lambda(j);
            end
        end
    end
else
    zd=zl;
    td=tl;
end



