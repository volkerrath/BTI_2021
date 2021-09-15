function [Rm,Rd] = ResMat(Jw,reg)
% calculates Generalized Inverse G
G=GI(Jw,reg);Rm=G*G';
if nargout>1, Rd=G'*G; end
end
