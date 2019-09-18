function [filelist]=mdir(mpath,mstrng)
% list[filelist]=mdir(mpath,mstrng) looks in directory mpath for files
% with specification mstrng. 
% Default for mstrng is '*", and for mpath is pwd.
% vr sep 2019
if nargin<2, mstrng='*'; end
if nargin<1, mpath=pwd; end
oldpath=pwd;
cd(mpath)
filelist=dir(mstrng);
cd(oldpath)