function [vn] = c2n(vc,d)
% VN=C2N(VC,D) interpolates cell center values  to nodes.
%
% function [vc] = c2n(vn,z) interpolates cellwise
% defined parameter to nodes. takes cell values vc and
% cell sizes d, gives nodal valuse vn. size of vn is
% length(vc)+1=length(d)+1.
%
% V. R.,  Aug.15, 2001 

nc=length(vc);nn=nc+1;
vn(2:nc) = (d(2:nc).*vc(2:nc)+d(1:nc-1).*vc(1:nc-1))./ ...
                     (d(2:nc)+d(1:nc-1));

vn(1) =vc(1) + d(1) *(vc(1)-vc(2))    /(d(1)+d(2));
vn(nn)=vc(nc)+ d(nc)*(vc(nc)-vc(nc-1))/(d(nc)+d(nc-1));
