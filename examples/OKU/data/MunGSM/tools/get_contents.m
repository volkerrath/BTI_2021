function [cont,dirflag] = get_contents(directory,varargin)
%  Get the contents of a specified directory
%
%   cont = get_contents(directory)
%   cont = get_contents(directory,options)
%   [cont,dirflag] = get_contents(...)
%
%   This function returns the contents of a specified
%   directory.
%
%   cont = get_contents(directory) returns the files and
%   folders in directory and returns them to the cell array
%   cont. It ignores hidden files and folders (those
%   starting '.'). Directory must be a character array
%   (string).
% 
%   cont = get_contents(directory,options) allow search
%   options to be specified. The options are specified in
%   parameter/value pairs, and include:
%       'rec'       {false} | true
%                   Search recursively within the
%                   subfolders of the specified directory.
%       'path'      {'relative'} | 'full'
%                   Specifies whether returned paths are
%                   full or relative to the specified
%                   directory.
%       'sort'      {false} | true
%                   Specify whether the output is sorted
%                   alphabetically.
%       'filter'    {'all'} | 'files' | 'folders' |
%                   '*.ext' | str
%                   This option allows a filter to be
%                   specified. 'files' returns names of all
%                   files in the directory. 'folders'
%                   returns names of all folders in the
%                   directory. '*.ext', where 'ext' is a
%                   user-specified file extension, returns
%                   all files with the extension '.ext'. str
%                   may be any string; only elements that
%                   contain str will be returned (files or
%                   folders). str is case-sensitive.
% 
%   [cont,dirflag] = get_contents(...) returns a logical
%   array dirflag, the same size as cont, indicating whether
%   each element is a directory.
% 
%   Examples
% 
%       Ex. 1
% 
%       % Return all m-files in the current directory
% 
%       cont = get_contents(cd,'filter','*.m')
% 
%       Ex. 2
% 
%       % Return all files in the current directory and its
%       % sub-directories
% 
%       cont = get_contents(cd,'rec',true)
% 
%       Ex. 3
% 
%       % Return all files in current directory with names
%       % containing 'foo'
%       
%       % may return files and folders:
%       [cont,dirflag] = get_contents(cd,'filter','foo')
% 
%       % use dirflag to limit:
%       cont = cont(~dirflag);

% !---
% ==========================================================
% Last changed:     $Date: 2013-06-01 10:14:12 +0100 (Sat, 01 Jun 2013) $
% Last committed:   $Revision: 248 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

% parse input arguments and arrange call(s) to 'main', which
% does the actual searching of directories

assert(ischar(directory),'directory must be a character array')

% Switch trap parses the varargin inputs
% default values
recflag = false;
pathflag = 'relative';
sortflag = false;
str = 'all';
% find values
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'path'
            pathflag=varargin{i+1};
        case 'rec'
            recflag=varargin{i+1};
        case 'sort'
            sortflag=varargin{i+1};
        case 'filter'
            str=varargin{i+1};
        otherwise
            error('Unknown option: %s\n',varargin{i});
    end
end

% check input options
assert(ischar(pathflag),'''path'' option must be a string')
assert(strcmp(pathflag,'relative') | strcmp(pathflag,'full'),...
    '''path'' option must ''relative'' or ''full''')
assert(islogical(recflag) & numel(recflag)==1,'''rec'' option must be logical')
assert(islogical(sortflag) & numel(sortflag)==1,'''sort'' option must be a logical')
assert(ischar(str),'str must be a character array')

% first pass: contents of top-level folder
[cont,dirflag] = main(directory,str);

% do the recursive bit, if recursion is requested
if recflag
    dirs = main(directory,'folders');
    count = length(dirs);
    n = 1;
    while n <= count % recursion requested
        [cont_temp,dirflag_temp] = main(dirs{n},str); % search them
        cont = [cont; cont_temp]; %#ok<AGROW> append search results
        dirflag = [dirflag; dirflag_temp]; %#ok<AGROW> append search results
        sdirs = main(dirs{n},'folders');
        dirs = [dirs; sdirs]; %#ok<AGROW>
        count = length(dirs);
        n = n+1;
    end
end

% remove full path
if strcmp(pathflag,'relative')
    if ~strcmp(directory(end),filesep)
        directory = [directory filesep];
    end
    for n = 1:length(cont)
        cont{n} = strrep(cont{n}, directory, '');
    end
end

% sort output (case insensitive)
if sortflag
    [~,IX] = sort(lower(cont));
    cont = cont(IX);
    dirflag = dirflag(IX);
end

% end of get_contents()


% ----------------------------------------------------------
% Local functions:
% ----------------------------------------------------------

% ----------------------------------------------------------
% main: get the contents
% ----------------------------------------------------------
function [cont,dirflag] = main(directory,str)

list = struct2cell(dir(directory));
dirbool = cell2mat(list(cellfun(@islogical,list(:,1)),:)); % return directory flags
list = list(1,:); % keep only file names
X = ~strncmp(list, '.', 1); % remove hidden files (those starting '.')
list = list(X);
list = list(:); % make column vector
dirbool = dirbool(X);
dirbool = dirbool(:); % make column vector

for n = 1:length(list)
    list{n} = fullfile(directory,list{n});
end

if nargin > 1
    % find filename extensions
    exts = cell(size(list));
    for n = 1:length(list)
        [~,~,exts{n}] = fileparts(list{n});
    end
    % filter
    if strncmp(str,'*.',2) % if extensions are requested
        ext = str(2:end);
        str = 'ext';
    end
    switch lower(str)
        case 'files'
            Y = ~dirbool;
        case 'folders'
            Y = dirbool;
        case 'ext'
            Y = strcmp(exts,ext);
        case 'all'
            Y = true(size(dirbool));
        otherwise % use literal search string
            Y = ~cellfun(@isempty,strfind(list,str));
    end
else
    Y = true(size(list));
end

% return search results
cont = list(Y);
dirflag = dirbool(Y);

% [EOF]
