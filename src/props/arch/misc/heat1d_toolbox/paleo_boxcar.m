function [T]=paleo_step(t,steptemp,steptime,debug)
% paleo_boxcar(t,stepamp,steptimedebug)
%
% initializes general step function for paleoclimate
% of amplitude stepamp at time steptime, given an input vector
% of temporal nodes at times t.
% if debug > 0, a control plot is produced.
%
% V. R., July 20, 2005
if nargin < 4, debug=0; end

nt=length(t);
T=steptemp(2)*ones(nt,1);   
maske=find(t<min(steptime));
T(maske)=steptemp(1);
maske=find(t>max(steptime));
T(maske)=steptemp(3);

if debug==1,
    figure;
    year2sec= 31557600;
    plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');
    set(gca,'XScale','log','XDir','reverse')
    grid on;ylim([-20 5]);xlim([10,100000]);
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    grid on;% legend('smoothed','original','Location','NorthWest')
    title(['test: set_paleo_boxcar'],'FontSize',14)
end

