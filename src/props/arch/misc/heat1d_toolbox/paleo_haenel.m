function [pt,it,T]=paleo_haenel(t,Tbase,debug)
% PALEO_HAENEL initializes paleoclimate after Haenel (1988)
%
% V. R.,  Sept. 25, 2001 
if nargin < 3, debug=0; end

year2sec=31557600;
th   = ...
- [100000, 75000, 70000, 65000, 55000, 35000, 26000, 15000, 10000, 9000, ...
   8000, 7200, 6700, 4500, 3500, 2480, 2000, 1580, 830, 680, 580, 380, ...
   330, 280, 230, 180, 130 ]*year2sec;
gth =  Tbase + ...
  [0, -4.9, -8., -4.9, -10.3, -6.2, -10.3, -7.0, 1.5, -2.0, 2.0, -1.0, 2.0, ...
  1.0, -0.5, 1.5, 1.5, 0.0,  0.8, 0.3, -0.5, -0.6, -0.7, -0.3, -0.6, -0.7, 0];

nth=length(th);
tt=t;it(tt<th(1))=1;
for j=1:nth-1
   lower=th(j);upper=th(j+1);
   step= tt>=lower & tt<upper;
   it(step)=j;
end
it(tt>=th(nth))=nth;
pt=gth;


if nargout > 3,
    T=gth(it);
    if debug==1,
        figure;
        plot(-t/year2sec',T, 'LineWidth',2,'Color','r');hold on;
        if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
        set(gca,'XScale','log','XDir','reverse')
        grid on;ylim([-20 5]);xlim([10,100000]);
        xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
        if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
        title(['test: set_paleo_haenel'],'FontSize',14)
    end
end
