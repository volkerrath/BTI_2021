%
%                             DISCLAIMER
%
%
% This file was generated by:
% ADiMat version 0.5.2 (*** DEBUG ***, beta, lcse, gcse(�), fwdmd, 2fwdmd, opt_ls, constfold, narg, vararg(�), Jun 27 2008) arch: i686-pc-linux-gnu
% compiled Jun 27 2008 with gcc 4.1.2 20070925 (Red Hat 4.1.2-33).
% Copyright 2001- 2007 Andre Vehreschild, Institute for
% Scientific Computing, Aachen University, D-52056 Aachen, Germany.
% http://www.sc.rwth-aachen.de/vehreschild/adimat/
% This file was augmented on Thu Oct 23 16:15:55 2008
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOLOOPSAVING -- Do not insert ls_* functions to encapsulate the loops over
%		 directional derivatives.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% VISITFWDMD  -- Use the visitor to generate the differentiated code.
% GRADFUNCPREFIX='sg_'

function [g_rhoi, rhoi]= sg_rhoiT(g_T, T)
   % calculate ice density in [kg/m**3]
   [n1, n2]= size(T); if n1== 1, g_tmp_rhoiT_00003= g_T' ;
      tmp_rhoiT_00003= T' ; % Update detected: T= some_expression(T,...)
      g_T= g_tmp_rhoiT_00003;
      T= tmp_rhoiT_00003;
   end
   % Identifier 'size' is ignored during differentiation.
   tmp_size_00000= size(T);
   g_tmp_size_00000= zeros(size(tmp_size_00000));
   g_rhoi= 0.00001* zeros(tmp_size_00000);
   rhoi= 0.00001* ones(tmp_size_00000); 
   clear tmp_size_00000 g_tmp_size_00000 ;
   tmp_rhoiT_00001= T< 0;
   g_tmp_rhoiT_00001= zeros(size(tmp_rhoiT_00001));
   g_rhoi(T< 0)= -0.151* g_T(tmp_rhoiT_00001);
   rhoi(T< 0)= 917.- 0.151* T(tmp_rhoiT_00001); 
   clear tmp_rhoiT_00001 g_tmp_rhoiT_00001 ;