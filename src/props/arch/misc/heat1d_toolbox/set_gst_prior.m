function [gst_apr,it]=set_gst_apriori(mode,gst_par)
% [gst_apr]=SET_GST_APRIORI(mode,gst_par) sets GST prior
% according to gstvale value of mode.

switch lower(mode)
   case {'constant'}
      n=gst_par.n;
      gt=gst_par.gst;
      gst_apr=gt*ones(1,n);
      it=ones(1,n);
   case {'read' 'file'}
      file=gst_par.file;t=gst_par.t;base=gst_par.base;
     [gstval step]=textread(file,'%f%f%*[^\n]',...
          'commentstyle','matlab','headerlines',header);
      ngstval=length(gstval);
      tt=t;it(tt<gstval(1))=1;
      for j=1:ngstval-1
          lower=gstval(j);upper=gstval(j+1);
          step= tt>=lower & tt<upper;
          it(step)=j;
      end
      it(tt>=gstval(ngstval))=ngstval;
      gst_apr=gstval;
   otherwise 
      error('set_gst_apriori: no mode given')
end   
