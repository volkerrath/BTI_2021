function [Rm] = ResMat(Jw,reg)
% calculates Generalized Inverse G
G=GenInv(Jw,reg);Rm=G*G';
if nargout>1, Rd=G'*G; end
end
