function ColorCodes = colorSpectrum(N)
% colorSpectrum Creates a spectrum of colors from red to green to blue. 
%   ColorCodes = colorSpectrum(N) ColorCodes contains N three-element RPG 
%   vectors. The RPG vectors follow the transition of red to green and 
%   green to blue. Avoids colors such as cyan and yellow, which are hard to
%   see in plots. N should be greater than 2.
% 
%   Ex. The following example shows the full spectrum made.
%       N = 1000;C=colorSpectrum(N);figure,hold on,
%       for i = 1:N,plot([1:10],ones(10,1)*i,'Color',C(i,:)),end,hold off

%   Copywrite 2013 Kirk T. Smith
%   

R = linspace(255,-255,N)';
G = [linspace(0,255,ceil(N/2)),linspace(255,0,ceil(N/2))]';
if mod(N,2)
    G = [G(1:ceil(N/2));G(ceil(N/2)+2:end)];
end
B = linspace(-255,255,N)';
ColorCodes = [R,G,B]/255;
ColorCodes = max(0,ColorCodes);
end


