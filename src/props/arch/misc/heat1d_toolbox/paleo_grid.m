function [pt,it]=paleo_grid(t,tbase,tstart,tend,ngrid)
% PALEO_GRID defines temporal inversion grid 
% (logarithmic discretization for paleoclimate)
%
% V.R., Dec 8, 2001
% year=31557600;
if nargin<4,ngrid=32;end
if nargin<3,tstart=abs(t(1));tend=abs(t(length(t)));end
if nargin<2,tbase=0;end


th=-logspace(log10(tstart),log10(tend),ngrid);
nth=length(th);
tt=t ;% -t(nt);
it(tt < th(1))=1;
for j=1:nth-1
      lower=th(j);  
      upper=th(j+1);
      step=tt>=lower & tt<upper;
      it(step)=j;      
end
it(tt >= th(nth))=nth;
pt=tbase*ones(1,nth);
