function [d,txt]=flexload(fname)
% A flexible text data file reader
% 
% Usage: [d,txt]=flexload(fname)
%
% This file code will read all the numbers found in a given file. It is
% flexible with whitespace, separators, etc. 
%
%
% INPUTS:
%   fname: filename of input file;
% 
% OUTPUTS:
%   d: a matrix with all the numbers in the file
%   txt: all the comments in the file. (prefixed with % or #);
%
% Aslak Grinsted 2014



fid=fopen(fname);
currows=0;
d=nan(10000,1);
txt=[];
while ~feof(fid)
    s=fgetl(fid);
    if isempty(s), continue, end
    ix=find((s=='%')|(s=='#')); %Comment Chars
    if ~isempty(ix)
        if nargout>1
            txt{end+1,1}=s(ix(1)+1:end);
        end
        s(ix(1):end)=[];
    end
    if isempty(s), continue, end

    v=regexpi(s,'([-+]?\d*\.?\d+(e[-+]?\d+)?|nan)','tokens'); %TODO... conjure up a nan handler 
    %v=cellfun(@(x)str2double(x{1}),v);  %slowest
    %v=str2double([v{1:end}]); %better
    v=cellfun(@(x)sscanf(x{1},'%f'),v); %better yet.
    
    currows=currows+1;
    
    nv=length(v);
    if nv>size(d,2)
        d(:,nv+1:end)=nan; %insert nans
    end
    if currows>size(d,1)
        d(end+(1:10000),:)=nan; %insert nans
    end

    d(currows,1:nv)=v;
    
end
fclose(fid);
d(currows+1:end,:)=[];

