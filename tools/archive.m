function zipfilename = archive(funcname,filename,archformat)
% archive - Creates a ZIP pr TGZ file containing all dependencies of a function
%
% zipfilename = archive(funcname)
% zipfilename = archive(funcname,filename)
% zipfilename = archive(funcname,filename, format)

% All M-files which are required by the specified function(s) and which are
% not inside MATLAB's "toolbox" directory.
% funcname is one or more function names, as a string or cell array of strings.
%
% If no filename is specified, the ZIP file is created in the current
% working directory with name X.zip where X is the first function name in
% the supplied list.  The file names in the archive file are relative to the
% common root directory directory of all the required files.  If no
% common root directory can be established (e.g. if files are on
% different drives) an error is thrown.
%
% Copyright 2006-2010 The MathWorks, Inc.

v=version;vnum=str2num(v(1)); 

if ~iscell(funcname)
    funcname = {funcname};
end

if isempty(funcname)
    error('No function names specified');
end

if ~iscellstr(funcname)
    error('Function names must be strings');
end


if nargin<2
    ffilename = fullfile(pwd,[ funcname{1} ]);
else 
    ffilename = fullfile(pwd,[ filename ]);
end

if nargin<3
    archformat='.zip';
end

archname=strcat([ffilename archformat]);

req = cell(size(funcname));
for i=1:numel(funcname)
    if vnum<8
    req{i} = mydepfun(funcname{i},1); % recursive
    else
    req{i} = matlab.codetools.requiredFilesAndProducts(funcname{i});
    end
    
end
req = vertcat(req{:}); % cell arrays of full file names
req = unique(req);

% Find the common root directory
d = i_root_directory(req);
% Calculate relative paths for all required files.
n = numel(d);
for i=1:numel(req)
    % This is the bit that can't be vectorised
    req{i} = req{i}(n+1:end); % step over last character (which is the separator)
end

switch lower(archformat)
    case {'.zip'}
      zip(archname,req,d);
    case {'tgz'}
      tar(archname,req,d);
    otherwise
      error(['archive: ',archformat,' unknown']);
end 
disp([' Created ',archname,' with ' num2str(numel(req)),' entries ']);


%%%%%%%%%%%%%%%%%%%%%
% Identifies the common root directory of all files in cell array "req"
function d = i_root_directory(req)

d = i_parent(req{1});
for i=1:numel(req)
    t = i_parent(req{i});
    if strncmp(t,d,numel(d))
        % req{i} is in directory d.  Next file.
        continue;
    end
    % req{i} is not in directory d.  Up one directory.
    count = 1;
    while true
        % Remove trailing separator before calling fileparts.  Add it
        % again afterwards.
        tempd = i_parent(d(1:end-1));
        if strcmp(d,tempd)
            % fileparts didn't find us a higher directory
            error('Failed to find common root directory for %s and %s',req{1},req{i});
        end
        d = tempd;
        if strncmp(t,d,numel(d))
            % req{i} is in directory d.  Next file.
            break;
        end
        % Safety measure for untested platform.
        count = count+1;
        if count>1000
            error('Bug in i_root_directory.');
        end
    end
end

%%%%%%%%%%%%%%%%%%%
function d = i_parent(d)
% Identifies the parent directory, including a trailing separator

% Include trailing separator in all comparisons so we don't assume that
% file C:\tempX\file.txt is in directory C:\temp
d = fileparts(d);
if d(end)~=filesep
    d = [d filesep];
end











