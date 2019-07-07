% Initialize MATLAB to use the opt_derivclass when computing ADiMat derivatives
% Part of the ADiMat tool: Automatic differentiation of Matlab programs
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% Aachen University of Technology, Andre Vehreschild, 2001, 2002
% Contact: vehreschild@sc.rwth-aachen.de
%
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

% Add the path of the ADiMat-toolbox  
adimat_home=getenv('ADIMAT_HOME');
if isempty(adimat_home)
   % Try to find the adimat executable and guess the location
   % of the derivclass. The methods are different for unix
   % and for Windows.
   if isunix
      % Look up adimat in the path
      [found,adimatpath]=unix('which adimat');
      % found is zero, if command was successful.
      if found==0
         % Extract the path from adimatpath. It is assumend, that
         % adimat resides in a directoy 'bin'.
         adimat_home=adimatpath(1:length(adimatpath)-11);
         if adimat_home(end)~='/'
            % The path is not valid. ADiMat resides in a directory
            % structure that is not GNU-conformant.
            adimat_home='';
         end
      else
         % The executable was not found. Assume default.
         % The default directory is '/usr/local/adimat-0.4-r9'.
         adimat_home='/usr/local/adimat-0.4-r9';
      end
      % Delete temporary variables.
      clear found adimatpath;
   else 
      if ispc 
         adimat_home=getenv('USERPROFILE');
         adimat_home=strcat(adimat_home, '/matlab/adimat/');
      else
         disp('Only unix is supported currently.');
      end
   end
end

if isunix
   tmp_full_path=strcat(adimat_home, '/share/adimat/opt_derivclass');
else
   % On windows the directory hierarchy is flatter.
   tmp_full_path=strcat(adimat_home, '/opt_derivclass');
end

if exist(tmp_full_path, 'dir')~=7
   % The directory does not exist. Ask the user for help.
         
   disp(['                     *** ADiMat- WARNING ***',10]);
   disp('The location of the support classes of ADiMat could not be determined.');
   disp('ADiMat augmented codes may not run. Set the environment varialbe');
   disp('"ADIMAT_HOME" to the location where the directory structure of ADiMat');
   disp(['can be found or add an addpath(',39,'your path to the derivclass',39,'); ']);
   disp('command to the Matlab startup files.');
else
   % The directory exists, set the path.
   addpath(tmp_full_path);
   setADoption(g_dummy, 'ADiMatHome', adimat_home);
   clear adimat_home;
end
clear tmp_full_path;

