function [rho]=rhofT(T,P,S)
% calculates water (brine) density depending on temperature (T, in C)
% pressure (P,in Pa), and salinity (S, in ppm) using the 
% formulae given by:
%  Mavko et al: Rock Physics Handbook, pp. 214ff, 1998.
% VR Feb. 3, 2003 
mode='brine'; 

if nargin < 3 mode='water'; S=0.; end
if nargin < 2 P=0.1; end
[n1,n2]=size(T); if n1~=1, T=T'; end 
[n1,n2]=size(P); if n1~=1, P=P'; end 
[n1,n2]=size(S); if n1~=1, S=S'; end 

% Pressure in MPa 
P=P/1.e6;P2=P.*P;
         T2=T.*T;T3=T2.*T;

switch lower(mode)
    case 'water'
        % pure water
        rhow=ones(1,n)+1.e-6* ...
            (-80*T-3.3*T2+0.00175*T3+489*P-2*T.*P+0.016*T2.*P-1.3e-5*T3.*P-(0.333-0.002*T).*P2);

        rho=rhow;
    case 'brine'
        % brine
        S=S/1.e6
        rhow=ones(1,n)+1.e-6* ...
            (-80*T-3.3*T2+0.00175*T3+489*P-2*T.*P+0.016*T2.*P-1.3e-5*T3.*P-(0.333-0.002*T).*P2);
        rho=rhow + S.* ...
            (0.668+0.44*S+1.e-6*(300*P-2400*P.*S+T.*(80+3*T-3300*S-13*P+47*P.*S)));
    otherwise
        rho=1000.;
end
rho=rho*1000;
