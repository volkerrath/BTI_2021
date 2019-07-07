function plot_m(par,z,xl,yl,labelx,labely,opts1,opts2)
% plot_m(par,d,xl,yl,labelx,labely) plots parameter + error bars 
% stored as columns of par as function of depth, given the intervals d.
% V.R.  Oct. 29, 2001 
if nargin < 7, opts1='b-';opts2='r--';end

stairs([par(1,1);par(:,1)],-z,opts1);hold on;grid on;
[n1,n2]=size(par);
if n2>1,
   for k=2:n2,
    stairs([par(1,k);par(:,k)],-z,opts2);hold on;
   end
end

if nargin>2, xlim(xl);ylim(yl);end
if nargin>4, xlabel(labelx);ylabel(labely);end
