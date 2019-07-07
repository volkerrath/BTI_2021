function [T]=paleo_osterkamp(t,Tbase,debug)
% PALEO_OSTERKAMP initializes paleotemperature of permfrost surface T in
% Prudhoe Bay Alaska.
%
% D. M.,  Mar. 20, 2003
if nargin < 3, debug=0; end
year2sec=31557600;
to   = ...
    - [240000, 210000, 200000, 185000, 150000, 135000, 130000, 120000, 90000, ...
    80000, 65000, 30000, 15000, 10000]*year2sec;
To =  Tbase + ...
    [-14.0, -12.0, -15.0, -13.0, -15.0,  -8.0, -11.0, -14.0, -12.0 -15.0, ...
    -13.0, -15.0, -8.0, -11.0];
nto=length(to);nt=length(t);
T(1:nt)=To(1);
for j=1:nto-1
    lower=to(j);upper=to(j+1);
    step= find(t>=lower & t<upper);
    T(step)=To(j);
end
last=find(t>=to(nto));
T(last)=To(nto);
if debug==1,
    figure;

    plot(-t/year2sec',T, 'LineWidth',2,'Color','r');hold on;
    if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
    set(gca,'XScale','log','XDir','reverse')
    grid on;ylim([-20 5]);xlim([10,100000]);
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
    title(['test: set_paleo_osterkamp'],'FontSize',14)
end
