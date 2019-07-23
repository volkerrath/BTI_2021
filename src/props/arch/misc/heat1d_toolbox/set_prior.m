function [mod]=set_prior(steptemp,steptime,L,tstart,tend,ngrid,debug)
% [mod]=set_prior(steptemp,steptime,L,tstart,tend,ngrid,debug)
%
% sets (posibly) smoothed prior model for paleoclimate inversion.
% for times > steptime temperature values steptemp arre assumed.
% the logarithmic inverse gride ranges from tstart to tend with
% ngrid steps. Sommotung is done by triangular fitering wit filter half
% length L.
% vr & dm, July 2005

year2sec=31557600;
steptime=[steptime 0];

th=-logspace(log10(tstart),log10(tend),ngrid);
nth=length(th);
mod(find(th < steptime(1)))=steptemp(1);
nm=length(mod);

if debug==1,

    disp(strcat(['POM (t<'  ,num2str(steptime(1))...
        ' ) Value: ',num2str(steptemp(1))]))
end

ns=length(steptime)-1;
for k=1:ns

    lower=steptime(k);
    upper=steptime(k+1);
    step=find(th>=lower & th<upper);
    mod(step)=steptemp(k);
    if debug==1,
        year2sec=31557600;
        disp(strcat(['Upper: ',num2str(upper/year2sec)...
            '   Lower: ',num2str(lower/year2sec)...
            '  Value: ',num2str(steptemp(k))]))
    end
end


mod(find(th >= steptime(ns)))=steptemp(ns);
nm=length(mod);
if debug==1,
    year2sec=31557600;
    disp(strcat(['Present (t> ',num2str(steptime(ns))...
    ' )  Value: ',num2str(steptemp(k))]))
end
%triangular window
if L ~= 0,
    w = tri(2*L+1); w = w/sum(w);
    X=mod';
    for k = L+1:nm-L
        aux = X(k-L:k+L);
        mod(k-L) = sum(w.*aux);
    end;
end
if debug==1,
    figure;
    [X,Y]=stairs(-th/year2sec,mod);
    plot(X,Y,'LineWidth',2,'Color','r','LineStyle','-');hold on;
    set(gca,'XScale','log','XDir','reverse')
    xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
    title(['test: set prior'],'FontSize',14)
    grid on;
end




