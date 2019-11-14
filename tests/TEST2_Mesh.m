function [ierr]=SITE_Mesh(name)
% Site-specific mesh generation

load common

ierr=0;
plotit= 0;
debug=0;

y2s=3600*24*365.25;s2y=1./y2s;

set_z = 1;
set_t = 1;

site            = 'TEST2';
zstart          = 0;
zend            = 5000;
z1max           = 2400;
nz              = 301;
dzstart=10; gfac=1.005;ngen=999;
ztype='log';

tstart= 110000*y2s;
tend= 30*y2s;
nt=401;
ttype= 'log';
dir= -1;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE TMP_MESH OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_Mesh_in.mat']);
if exist(F)
    disp([' ...' mfilename ' defaults overwritten from ', F])
    load(F); mstruct(mesh_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if plotit
    set_graphpars
end

if set_z
    disp(['   ']);
    disp(strcat([ ' ...Set up meshes for site ', site]));
    % % SPATIAL MESH
    switch lower(ztype)
        case{'special' 'mixed'}
            disp([' ...set up ' ztype ' temporal mesh ']);
            z1=linspace(zstart,z1max,nz);
            dzstart=8.; gfac=1.02;ngen=500;dzn=dzstart*gfac.^[1:ngen];
            z2=cumsum(dzn);
            z=[z1 z1max+z2]; z=z(z<zend); z=[z zend];
        case{'read' 'inp'}
            meshfile='Input_zmesh';
            disp([' ']);disp([' ...load  spatial mesh from ',meshfile]);
            load(meshfile);
            z=z(z<zend); z=[z zend];
        case{'linear' 'lin'}
            z=linspace(zstart,zend,nz);
        case {'logarithmic','log'}
            disp([' ...set up ' ztype ' spatial mesh ']);
%             dzn=dzstart*gfac.^[1:ngen];
            z=logspace(log10(zstart),log10(zend),nz-1);
            z=[0 z]; 
        otherwise
            disp([' ']);disp([' ...spatial mesh ',ztype,' not implemented!']);
    end
    z = unique(z,'stable');
    dz=diff(z); nz=length(z);
    ip=[1:nz-1];zm=0.5*(z(1:nz-1)+z(2:nz));
    
    F=strcat([name '_DepthGrid.mat']);
    disp([' ...spatial mesh written to: ', F]);
    save(F,'z','dz','nz','ip','zm')
end

if set_t
    % TEMPORAL MESH
    switch lower(ttype)
        case{'special' 'mixed'}
            disp([' ...set up ' ttype ' temporal mesh ']);
            t1=[-110000:100:0]*y2s;
            t2=[-20000:50:-3000]*y2s;
            t3=[-3000:20:0]*y2s;
            t=union(t1(:),union(t2(:),t3(:)));t=sort(t);
        case{'read'}
            meshfile='Input_tmesh';
            disp([' ...load  temporal mesh from ',meshfile]);
            load(meshfile);
        case {'logarithmic','log'}
            disp([' ...set up ' ttype ' temporal mesh ']);
            [t,dt]= set_mesh(tstart, tend, nt, ttype, dir, debug);
        case{'linear' 'lin'}
            disp([' ...set up ' ttype ' temporal mesh ']);
            t=linspace(tstart,tend,nt);

        otherwise
            disp([' ']);disp([' ...remporal mesh ',ttype,' not implemented!']);
    end
    t = unique(t,'stable');
    dt=diff(t);nt=length(t);it=[1:nt];
    tm=0.5*(t(1:nt-1)+t(2:nt));
    
    F=strcat([name '_TimeGrid.mat']);
    save(F,'t','dt','nt','it','tm')
    disp([' ...temporal mesh written to: ', F]);
    %
    disp([mfilename '    spatial mesh: ',num2str(nz),' temporal mesh:',num2str(nt)]);
end
