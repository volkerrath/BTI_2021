function plot_d(d,z,xl,yl,labelx,labely,opts1,opts2)
% plot_d(d,dz,xl,yl,labelx,labely) plots (multiple) data
% stored as columns of v as function of depth, given the intervals d.
% v.r.  Aug.19,2001 
if nargin < 7, opts1=['b-'];opts2=['r--'];end

plot(d(:,1),-z,opts1);hold on;grid on;
[n1,n2]=size(d);
if n2>1,
   for k=2:n2,
      plot(d(:,k),-z,opts2);hold on;
   end
end

if nargin>2,xlim(xl);ylim(yl);end
if nargin>4, xlabel(labelx);ylabel(labely);end
