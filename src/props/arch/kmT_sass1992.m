function kmT = kmT_sass1992(k,T,p)
% SASS1992 Temperature dependent thermal conductivity
%
% L = SASS1992(p,T,L25)
%
% L25: lambda_25, thermal conductivity at room temperature
% T: vector of temperature values
% p: 1x3 vector of parameters a,b,c
%
% L = L0 / ( a + T * (b - c/L0) )

% aha, 5 jan 2005
% kb based on Vosteen (2003)
T=T(:);k=k(:);
kb = 0.54*k + 0.5 * sqrt( 1.16*k.^2 - 0.39*k );
kmT = kb ./ ( p(:,1) + T .* ( p(:,2) - p(:,3)./kb) );


