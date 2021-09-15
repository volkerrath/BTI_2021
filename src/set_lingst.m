function [S,T]=set_lingst(t,amp,tim,L,debug)
% [S,T]=SET_GST(TIME,STEPAMP,STEPTIME,L,DEBUG) initializes general step function
% for paleoclimate with amplitude stepamp at time steptime, g
% iven an input vector of temporal nodes at times t. L is the 
% length of a smoothing operator applied to the temperature time series.
% if debug > 0, a control plot is produced.
%
% V. R., July 20, 2005
if nargin < 5, debug=0; end;
if nargin <4, L=0; end;
t=t(:); [nt,n0]=size(t);


method = 'linear';

meanamp=mean(amp);

T = interp1(tim,amp,t,method);
ntf=find(isfinite(T),1,'first');T(1:ntf-1)=T(ntf);
ntl=find(isfinite(T),1,'last'); T(ntl+1:length(T))=T(ntl);


S = T(:);


if L ~= 0,
    T=T';
    nT=length(T);
    % triangular window
    w = tri(2*L+1); w = w/sum(w);w=w';

    for k = L+1:nT-L
        aux = T(k-L:k+L);
        S(k) = sum(w.*aux);
    end;
end

if debug >0,
    figure;
    year2sec= 31557600;
    plot(-t/year2sec',S, 'LineWidth',2,'Color','r');hold on;
    if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
    set(gca,'XScale','log','XDir','reverse')
    grid on;
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
    title(['test: set_paleo_boxcar'],'FontSize',14)
    if debug >=1,
        fielname='GST.ps';
        saveas(gcf,filename,'psc2')
    end
end






