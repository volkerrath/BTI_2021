function [pt,it,c,l,u]=set_mgsth(t,tbase,tstart,tend,ngrid,gmode,thin)
% defines temporal inversion grid
% (discretization for paleoclimate)
%
% V.R., Nov, 2014
% V.R., June, 2013
if nargin<6,gmode='log';end
if nargin<5,ngrid=32;end
if nargin<3,tstart=abs(t(1));tend=abs(t(length(t)));end
if nargin<2,tbase=0;end
if nargout > 2
    c=NaN*ones(1,ngrid-1);u=c;l=c;
end


switch lower(gmode)
    case {'log' 'logarithmic'}
        th=-logspace(log10(tstart),log10(tend),ngrid);
        nth=length(th);
        tt=t ;% -t(nt);
        it(tt < th(1))=1;
        for jj=1:nth-1
            low=th(jj);
            upp=th(jj+1);
            step=tt>=low & tt<upp;
            it(step)=jj;
            if nargout > 2
                if ~isempty(step),
                    c(jj)=mean(log(tt(step)));
                    l(jj)=tt(find(step,1,'first'));
                    u(jj)=tt(find(step,1,'last'));
                end
            end
        end
        it(tt >= th(nth))=nth;
        pt=tbase*ones(1,nth);
        if nargout > 2, c=exp(c); end
    case {'lin','linear'}
        th=-linspace(tstart,tend,ngrid);
        nth=length(th);
        tt=t ;% -t(nt);
        it(tt < th(1))=1;
        for jj=1:nth-1
            low=th(jj);
            upp=th(jj+1);
            step=tt>=low & tt<upp;
            it(step)=jj;
            if nargout > 2
                if ~isempty(step),
                    c(jj)=mean(tt(step));
                    l(jj)=tt(find(step,1,'first'));
                    u(jj)=tt(find(step,1,'last'));
                end
            end
        end
        it(tt >= th(nth))=nth;
        pt=tbase*ones(1,nth);
    case {'in','input'}
        th=thin;
        nth=length(th);
        tt=t ;% -t(nt);
        it(tt < th(1))=1;
        for jj=1:nth-1
            low=th(jj);
            upp=th(jj+1);
            step=tt>=low & tt<upp;
            it(step)=jj;
            if nargout > 2
                if ~isempty(step),
                    c(jj)=mean(tt(step));
                    l(jj)=tt(find(step,1,'first'));
                    u(jj)=tt(find(step,1,'last'));
                end
            end
        end
        it(tt >= th(nth))=nth;
        pt=tbase*ones(1,nth);    
   otherwise
        disp('Unknown method.')
end
