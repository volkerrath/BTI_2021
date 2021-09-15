function res= option(cmd, val)
% ADDERIV/PRIVATE/OPTION -- Set and get options.
%
% Copyright 2003, 2004 Andre Vehreschild, Inst. f. Scientific Computing
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

persistent options

if isempty(options)
   options= clearAll(options);
end

if nargin>1
   % Set data !!!
   switch cmd
   case 'NumberOfDirectionalDerivatives'
      if options.NumDirDeriv==0
         options.NumDirDeriv= val;
      elseif options.NumDirDeriv~=val
         error('All derivative object have to store the same number of directional derivatives.');
      end        
   case 'ADiMatHome'
      options.ADiMatHome= val;
   case 'Version'
       error('The version of ADiMat can not be modified.');
   otherwise
      error('Unknown option!');
   end
else
   % Get data !!!
   switch cmd
   case 'NumberOfDirectionalDerivatives'
      res= options.NumDirDeriv;
   case 'SubsRefBug'
      res= options.subsrefbug;
   case 'ADiMatHome'
      res= options.ADiMatHome;
   case 'Version'
      res= options.Version;
   case 'ClearAll'
      options.NumDirDeriv= 0;
      warning('ADiMat:clearOptionswarning', ...
              'All AD-options are set back to defaults.');
      res= [];
   case 'DerivativeClassName'
      res='adderiv';
   case 'DerivativeClassVersion'
      res=0.5;
   case 'DerivativeClassKind'
      res='adderiv 0.5';
   otherwise
      error('Unknown option!');
   end
end

function options= clearAll(options)
   adimat_version= 0.4 ;
   adimathome= getenv('ADIMAT_HOME');  
   if ispc & isempty(adimathome)
       try
          adimathome= winqueryreg('HKEY_LOCAL_MACHINE', ['SOFTWARE\ADiMat\' num2str(version)], 'home');
       catch 
          adimathome= [];
       end
       if isempty(adimathome)
           adimathome= which('adimat.exe');
           if ~ isempty(adimathome)
              adimathome= adimathome(1:end-10);    
           else
              warning('ADiMat:HomeWarning', 'Can not find home of ADiMat. ADiMat is not properly installed.');
           end
       end
   end
   
   mlrel=str2double(version('-release'));
   if isnan(mlrel)
      % Mathworks unfortunately changed the release numbering scheme
      % starting with MATLAB 2006a. Fortunately the bug is not present there.
      srefbug= false;
   else
      srefbug=mlrel<14.0;
   end
   options=struct('ADiMatHome', adimathome, 'Version', adimat_version, 'NumDirDeriv', 0, 'subsrefbug', srefbug);

% vim:sts=3:
