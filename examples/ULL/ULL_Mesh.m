function [ierr]=SITE_Mesh(name)
% Site-specific mesh generation

load common

ierr=0;
plotit= 0;
debug=1;

y2s=3600*24*365.25;s2y=1./y2s;

site            = 'ULL';
zstart          = 0;
zend            = 5000;
zlmax           = 1540;
nz              = 301;
dzstart=4; gfac=1.005;ngen=999;
ztype='log';


tstart= 110000*y2s;
tend= 10*y2s;
nt=357;
ttype= 'log';
dir= -1;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE OKU_MESH OVERWRITE DEFAULTS ABOVE!h
F=strcat([name,'_Mesh_in.mat']);
if exist(F)
    load(F); mstruct(mesh_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if plotit
    set_graphpars
end

disp(['   ']);
disp(strcat([ ' ...Set up meshes for site ', site]));
% % SPATIAL MESH
switch lower(ztype)
    case{'special' 'mixed'}
        z1=linspace(zstart,zlmax,nz);
        dzstart=8.; gfac=1.02;ngen=500;dzn=dzstart*gfac.^[1:ngen];
        z2=cumsum(dzn);
        z=[z1 zlmax+z2]; z=z(z<zend); z=[z zend];
        dz=diff(z);nz=length(z);
        ip=[1:nz-1]; zm=0.5*(z(1:nz-1)+z(2:nz));
    case{'read' 'inp'}
        meshfile='Input_zmesh';
        disp([' ']);disp([' ...load  spatial mesh from ',meshfile]);
        load(meshfile);
        z=z(z<zend); z=[z zend];
        dz=diff(z);nz=length(z);
        ip=[1:nz-1]; zm=0.5*(z(1:nz-1)+z(2:nz));
    case {'log'}
        disp([' ']);disp([' ...set up ' ztype ' spatial mesh ']);
        dzn=dzstart*gfac.^[1:ngen];
        z=[0 cumsum(dzn)]; z=z(z<zend); z=[z zend];
        dz=diff(z); nz=length(z);
        ip=[1:nz-1]; zm=0.5*(z(1:nz-1)+z(2:nz));
    otherwise
        disp([' ']);disp([' ...spatial mesh ',ztype,' not implemented!']);
end

F=strcat([name '_DepthGrid']);
disp([' ...spatial mesh written to: ', F]);
save(F,'z','dz','nz','ip','zm')



% T% % TEMPORAL MESH
switch lower(ttype)
    case{'special' 'mixed'}
        t1=[-110000:100:0]*y2s; 
        t2=[-20000:50:-3000]*y2s;   
        t3=[-3000:20:0]*y2s;
        t=union(t1(:),union(t2(:),t3(:)));t=sort(t);it=[1:nt];
        dt=diff(t);nt=length(t);
    case{'read'}
        meshfile='Input_tmesh';
        disp([' ']);disp([' ...load  temporalal mesh from ',meshfile]);
        load(meshfile);
        dt=diff(t);nt=length(t);it=[1:nt];
    case {'log'}
        [t,dt]= set_mesh(tstart, tend, nt, ttype, dir, debug);
        dt=diff(t);nt=length(t);it=[1:nt];
        tm=0.5*(t(1:nt-1)+t(2:nt));
    otherwise
        disp([' ']);disp([' ...spatial mesh ',ztype,' not implemented!']);
end

F=strcat([name '_TimeGrid']);
save(F,'t','dt','nt','it','tm')
disp([' ...temporal mesh written to: ', F]);

end
