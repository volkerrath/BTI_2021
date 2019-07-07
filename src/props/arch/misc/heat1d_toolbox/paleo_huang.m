function [pt,it,T]=paleo_huang(t,Tbase,debug)
% PALEO_HUANG initializes global warming after Huang & Pollack (2001) 
%
% V. R.,  Oct. 10, 2001 
if nargin < 3, debug=0; end
year=31557600;
tp  = ...
   - [500 400 300 200 100 0.0001];  
gtp = Tbase + ...
   [ -.994 -.959   -.886 -.751 -.525  0;  ...
    -1.098 -1.046  -.959 -.810 -.566  0;
     -.890 -.872   -.813 -.692 -.484  0];
 nt=length(t);nth=length(tp);
tt=t;it(tt<tp(1,1))=1;
for j=1:nth-1
   lower=thp(j);upper=tp(j+1,1);
   step= tt>=lower & tt<upper;
   it(step)=j;
end
it(tt>=tp(nth))=nth;
pt=gtp(:,1);

if nargout > 3,
    T=gtp(it);
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
end
