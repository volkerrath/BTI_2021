function [S]=set_pntgst(t,tpoint,apoint,method,debug)
% [S,T]=SET_GST(TIME,TPOINT,APOINT,METHOD,DEBUG) initializes gsth from
% step function for paleoclimate with amplitude stepamp at time steptime, 
% given an input vector of temporal nodes at times t. L is the 
% length of a smoothing operator applied to the temperature time series.
% if debug > 0, a control plot is produced.
% Last change: vr, July 20, 2019
if nargin < 4, debug=0; end
logint = 1;

if logint
    tpoint = log(abs(tpoint)) ;
    t=log(abs(t));
end

S=interp1(tpoint,apoint,t,method);
% S(t<tpoint(1))=S(1);
% S(t>tpoint(length(tpoint)))=S(length(tpoint));

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
        filename='GST.png';
        saveas(gcf,filename,'png')
    end
end






