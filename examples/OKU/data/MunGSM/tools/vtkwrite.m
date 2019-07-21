function vtkwrite( filename,dataType,varargin )
% VTKWRITE Write 3D Matlab array into VTK file format.
%  vtkwrite(filename,'STRUCTURED_GRID',x,y,z,'vector',title,u,v,w) writes a
%  3D vector field data into VTK file whose name is specified by the string
%  filename. (u,v,w) are the vector components at the points (x,y,z). The
%  matrices x,y,z,u,v,w must all be the same size and contain
%  corrresponding position and vector component. The string title specifies
%  the name of the vector field to be saved. 'STRUCTURED_GRID' specify the
%  type of the dataset to be saved. Currently it's the only type supported
%  by the function. While redundant at this moment, it is written this way
%  for the ease of future expansion of support for other types. 
%
%  vtkwrite(filename,'STRUCTURED_GRID',x,y,z,'scalar',title,r) writes a 3D
%  scalar field data into VTK file whose name is specified by the string
%  filename. r is the scalar value at the points (x,y,z). The matrices
%  x,y,z,r must all be the same size and contain the corresponding position
%  and scalar values. 
%
%  vtkwrite(filename,'STRUCTURED_GRID,x,y,z,'vector',title1,u,v,w,'scalar',
%  title2,r) writes a 3D field that contains both vector and scalar values.
%  x,y,z,u,v,w must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%
%  Version 1.0.0
%  Copyright 2014 Chaoyuan Yeh
%

fid = fopen(filename, 'w'); 
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'VTK from Matlab\n');
switch upper(dataType)
    case 'STRUCTURED_GRID'
        fprintf(fid, 'BINARY\n');
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        if numel(varargin)<6, error('Not enough input arguments'); end
        if sum(size(x)==size(y) & size(y)==size(z))~=3
            error('Input dimesions do not match')
        end
        n_elements = numel(x);
        fprintf(fid, 'DATASET STRUCTURED_GRID\n');
        fprintf(fid, ['DIMENSIONS ' num2str(size(x,2)) ' ' num2str(size(x,1)) ' ' num2str(size(x,3)) '\n']);
        fprintf(fid, ['POINTS ' num2str(n_elements) ' float\n']);
        fwrite(fid, [x(:)';y(:)';z(:)'],'float','b');
        fprintf(fid, ['\nPOINT_DATA ' num2str(n_elements)]);
        
        % Parse remaining argument.
        vidx = find(strcmpi(varargin,'vector'));
        sidx = find(strcmpi(varargin,'scalar'));
        if vidx~=0
            for ii = 1:length(vidx)
                title = varargin{vidx(ii)+1};
                fprintf(fid, ['\nVECTORS ', title,' float\n']);
                fwrite(fid, [reshape(varargin{vidx(ii)+2},1,n_elements);...
                reshape(varargin{vidx(ii)+3},1,n_elements);...
                reshape(varargin{vidx(ii)+4},1,n_elements)],'float','b');
            end
        end
        if sidx~=0
            for ii = 1:length(sidx)
                title = varargin{sidx(ii)+1};
                fprintf(fid, ['\nSCALARS ', title,' float\n']);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                fwrite (fid, reshape(varargin{vidx(ii)+2},1,n_elements),'float','b');
            end
        end
end
fclose(fid);        
end