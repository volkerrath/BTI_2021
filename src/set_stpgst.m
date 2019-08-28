function [S,T]=set_stpgst(t,stepamp,steptime,L,pom,debug)
% [S,T]=SET_GST(TIME,STEPAMP,STEPTIME,L,DEBUG) initializes gsth from
% step function for paleoclimate with amplitude stepamp at time steptime, 
% given an input vector of temporal nodes at times t. L is the 
% length of a smoothing operator applied to the temperature time series.
% if debug > 0, a control plot is produced
% Last change: vr July 20, 2019
if nargin < 6, debug=0; end
if nargin < 5, pom=NaN; end
if nargin < 4, L=0; end

ns=length(steptime); nt=length(t);
if ~isfinite(pom)
    pom=stepamp(ns)-4;
end

maske=find(t<steptime(1));
% maske
T(maske)=pom;
for i=2:ns 
%     i
    maske=find(t<steptime(i) & t>=steptime(i-1));
%     maske
    T(maske)=stepamp(i-1);
end
maske=find(t>steptime(ns));
T(maske)=stepamp(ns);
S=T(:);

if L ~= 0
    T=T';
    nT=length(T);
    % triangular window
    w = tri(2*L+1); w = w/sum(w);

    for k = L+1:nT-L
        aux = T(k-L:k+L);
        S(k) = sum(w.*aux);
    end
end

if debug >0
    figure;
    year2sec= 31557600;
    plot(-t/year2sec',S, 'LineWidth',2,'Color','r');hold on;
    if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
    set(gca,'XScale','log','XDir','reverse')
    grid on;
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
    title(['test: set_paleo_boxcar'],'FontSize',14)
    if debug >=1
        fielname='GST.png';
        saveas(gcf,filename,'png')
    end
end






