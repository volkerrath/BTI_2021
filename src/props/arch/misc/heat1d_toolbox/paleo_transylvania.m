function [pt,it,T]=paleo_transylvania(t,tbase,debug)

year2sec= 31557600;
gt_prior=[10.0   10.0    10.0     7.5     7.5 ...
    10.5    10.5     9.75    9.75    11.25 ...
    11.25    9.75    9.75   -2.25    -2.25 ...
    0.75    0.75    9.75]'+tbase;
t_prior=-[ 0.0    50.0   100.0   150.0    500.0 ...
    700.0  1000.0  1050.0  2000.0   2050.0 ...
    7000.0  7050.0 10000.0 10100.0  35000.0 ...
    35100.0 65000.0 65100.0]'*year2sec;
t_prior=flipud(t_prior);gt_prior=flipud(gt_prior);
nt_prior=length(t_prior);
tt=t;it(tt<t_prior(1))=1;
for j=1:nt_prior-1
    lower=t_prior(j);upper=t_prior(j+1);
    step= tt>=lower & tt<upper;
    it(step)=j;
end
it(tt>=t_prior(nt_prior))=nt_prior;
pt=gt_prior;

if nargout > 3,
    T=gt_prior(it);
    if debug==1,
        figure;

        plot(-t/year2sec',T, 'LineWidth',2,'Color','r');hold on;
        if L ~= 0, plot(-t/year2sec',T, 'LineWidth',2,'Color','b','LineStyle','--');end
        set(gca,'XScale','log','XDir','reverse')
        grid on;ylim([-20 5]);xlim([10,100000]);
        xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
        if L ~= 0, grid on;legend('smoothed','original','Location','NorthWest');end
        title(['test: set_paleo_transylvania'],'FontSize',14)
    end
end
