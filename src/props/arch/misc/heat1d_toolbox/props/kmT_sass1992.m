function l = kmT_sass1992(p,T,l25)
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

% l0 based on Vosteen (2003)
l0 = 0.54*l25 + 0.5 * sqrt( 1.16*l25.^2 - 0.39*l25 );

l = l0 ./ ( p(:,1) + T .* ( p(:,2) - p(:,3)./l0) );


