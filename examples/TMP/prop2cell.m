function Ks=prop2cell(param,zparam,z,Top,Bottom,avgmeth)
% BULK THERMAL CONDUCTIVITY
if nargin < 6, avgmeth = 'h';end

z = z(:);
zm = 0.5*(z(1:end-1)+z(2:end));
nc = length(zm);
Ks=NaN(size(zm));

for cell = 1:nc
    incell = find(zparam >= z(cell) & zparam<z(cell+1));
    if ~isempty(incell)
        Kc=param(incell);w=ones(size(Kc));
        Ks(cell) = vavg(Kc,w,avgmeth);
    else
        %         disp(strcat(['...>>> cell: ',num2str(cell),' centered at ',...
        %                      num2str(zm(cell)),' m empty']))
    end
end

uppval=find(isfinite(Ks),1,'first');
if ~exist('Top','var'), Top=Ks(uppval);end
Ks(1:uppval-1)=Top;

lowval=find(isfinite(Ks),1,'last');
if ~exist('Bottom','var'),Bottom=Ks(lowval);end
Ks(lowval+1:nc)=Bottom;

okval=find(isfinite(Ks));zi=zm(okval);Ki=Ks(okval);
Ks=interp1(zi,Ki,zm);