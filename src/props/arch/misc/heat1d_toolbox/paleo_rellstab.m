function [Ts]=paleo_haenel(t,Tbase,debug)
% PALEO_HAENEL initializes paleoclimate after Haenel (1988)
%
% V. R.,  Sept. 25, 2001
if nargin < 3, debug=0; end
year=31557600;
th   = ...
    - [1000000, 1000000, 900000, 800000, 600000, 450000, 300000, 250000, 200000, 125000, 70000, ...
    40000, 25000, 14000, 10000, 8000, 5000, 3200, 2000, 1000, 850, 650, 400, ...
    280, 230, 200, 170, 80, 50]*year;
gth =  Tbase + ...
    [0, -10.0, -5.0, -10.0, -7.0, -2.0, -13.0, -2.0, -13.0, -2.0, -11.0, -4.0, -13.0, -6.0, 0.0, ...
    1.5, 0.5, -0.5, 0.0, -0.3,  0.4, -0.7, -1.0, -0.8, -0.7, -0.6, -0.8, -0.4, 0];

nth=length(th);
tt=t;it(tt<th(1))=1;
for j=1:nth-1
    lower=th(j);upper=th(j+1);
    step= tt>=lower & tt<upper;
    it(step)=j;
end
it(tt>=th(nth))=nth;
Ts=gth(it);

if debug==1,
    figure;

    plot(-t/year2sec',T, 'LineWidth',2,'Color','r');hold on;
    if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
    set(gca,'XScale','log','XDir','reverse')
    grid on;ylim([-20 5]);xlim([10,100000]);
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
    title(['test: set_palew_huang'],'FontSize',14)
end
