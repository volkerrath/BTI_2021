function [err]=put_data(file,data,text)
% PUT_DATA writes DATA and ERRORS 
% values to file
%
% Arguments:
% file                  filename for temperature log
% text                  header             
% data                  structure associated with data       
%
%
%% V. R.,Nov. 3, 2002 

if nargin<3, text=strcat([file ' written on ' date]); end

err=0;
header=0;
fid=fopen(file,'w');
fprintf(fid,'%% %s\n',text);

z=data.z;T=data.T;Terr=data.Err;


fprintf(fid,'%12g %12g %12g \n',...     
         [z;T;Terr]);

                                                                   
fclose(fid);
