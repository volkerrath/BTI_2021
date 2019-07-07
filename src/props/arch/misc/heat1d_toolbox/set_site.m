function [data,model]=set_site(site,z)

% prepare data for site - special version OUTUKUMPU
tit ='OKU';
if site ~= tit,
    message=strcat([ 'Site ',site,' not equal ',tit,'. Stopped.'])
    error(message)
end
disp(strcat([ 'Preprocessing site ', tit]));
step=0;

nz=length(z);

% READ MEASUREMENTS AND DATA FOR BOREHOLE

% GENERAL SETTINGS
out   = 0;
debug = 0;
graph = 0;
estq  = 0;
% TEMPERATURE ERROR
err=0.1;
% DZ AVERAGE TEMPERATURE
dzT=3;
% DZ AVERAGE THERMAL CONDUCTIVITY
dzL=10;
% FILTER
windowSize = 3;
% PHYSICAL PARAMETERS 
por = 0.01;
GST0  = 5.2;
Q0  =42.e-3;

if estq, zEst =[2000 2100 2200 2300 2400]; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP I: READ & PREPROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ 'Step ',num2str(step),': read data - smooth and average']));


OT    =   importdata('OKU_Temp_orig.dat');
OL    =   importdata('OKU_Lamb_orig.dat');

zTo = OT(:,1);
To = OT(:,2);
zLo = OL(:,1);
Lo = OL(:,2);



zgT=0:dzT:max([max(zTo),max(zLo)]);
% TEMPERATURE
for cell = 1:length(zgT)-1
    zTog(cell) = 0.5*(zgT(cell)+zgT(cell+1));
    incell = find(zTo >= zgT(cell) & zTo<zgT(cell+1));
    Toa(cell) = mean(To(incell));
    Tom(cell) = median(To(incell));
    Tot(cell)=  mad(To(incell,1),1);
    Tov(cell)=  var(To(incell,1),1);
    incT(cell) = length(incell);
end
% THERMAL CONDUCTIVITY
zgL=0:dzL:max([max(zTo),max(zLo)]);
for cell = 1:length(zgL)-1
    zLog(cell) = 0.5*(zgL(cell)+zgL(cell+1));
    incell = find(zLo >= zgL(cell) & zLo<zgL(cell+1));
    Loa(cell) = mean(Lo(incell));
    Lom(cell) = median(Lo(incell));
    Lot(cell)=  mad(Lo(incell,1),1);
    Lov(cell)=  var(Lo(incell,1),1);
    incL(cell) = length(incell);
end

% FILTER

Tomf=filter(ones(1,windowSize)/windowSize,1,Toa);
Totf=filter(ones(1,windowSize)/windowSize,1,Tot);
Lomf=filter(ones(1,windowSize)/windowSize,1,Lom);
Lotf=filter(ones(1,windowSize)/windowSize,1,Lot);

definx=isfinite(Tomf);zx=zTog(definx);Tx=Tomf(definx);
GradT=diff(Tx)./diff(zx);Gradz= 0.5*(zx(1:length(zx)-1)+zx(2:length(zx)));
definx=isfinite(Lomf);zx=zLog(definx);Lx=Lomf(definx);
LInt = interp1(zx,Lx,Gradz);
Qobs = LInt.*GradT;

if graph,
    % PLOT TEMPERATURES & GRADIENTS
    figure
    plot(To,zTo,'-b','LineWidth',1); hold on
    plot(Tomf,zTog,'-r','LineWidth',1); hold on
    % plot(Tomf+Totf,zTog,'or','LineWidth',2); hold on
    % plot(Tomf-Totf,zTog,'or','LineWidth',2); hold on
    xlim([0 50]);
    ylim([0 2600]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWehttp://legacy.orie.cornell.edu/~davidr/ight',fontwg)
    title(strcat(['Outukumpu - Sept 5, 2008 (filtered,median \pm mad)']),'FontSize',fontsz,'FontWeight',fontwg);
    S1=strcat(['original data']);
    legend(S1,'location', 'northeast')
    grid on
    file=strcat(['OKU1_temp.eps']);
    saveas(gcf,file,'epsc2')

    figure
    %plot(To,zTo,'-b','LineWidth',1); hold on
    plot(GradT*1000,zGrad,'-r','LineWidth',2); hold on
    % plot(Tomf+Totf,zTog,'or','LineWidth',2); hold on
    % plot(Tomf-Totf,zTog,'or','LineWidth',2); hold on
    %xlim([0 50]);
    ylim([0 2600]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('\nabla T (K/km)','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['Outukumpu - Sept 5, 2008 (filtered,median \pm mad)']),'FontSize',fontsz,'FontWeight',fontwg);
    S1=strcat(['original data']);
    legend(S1,'location', 'southwest')
    grid on
    file=strdisp(strcat([ 'Step 2: interpolation to spatial grid']));cat(['OKU1_tgrd.eps']);
    saveas(gcf,file,'epsc2')

    % PLOT CONDUCTIVITIY
    figure
    plot(Lo,zLo,'-b','LineWidth',1); hold on
    plot(Lomf,zLog,'-r','LineWidth',2); hold on
    plot(Lomf-Lotf,zLog,':g','LineWidth',2); hold on
    plot(Lomf+Lotf,zLog,':g','LineWidth',2); hold on

    ylim([0 2600]);
    xlim([0 10]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('\lambda (W m^{-1}K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['Outukumpu - Sept 5, 2008 (filtered median \pm mad)']),'FontSize',fontsz,'FontWeight',fontwg)
    S1=strcat(['original data']);
    S2=strcat(['median, filtered']);
    legend(S1,S2,'location', 'northeast')
    grid on
    file=strcat(['OKU1_Lamb.eps']);
    saveas(gcf,file,'epsc2')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LInt*GradT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP II: INTERPOLATE TO GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step=step+1;
disp(strcat([ 'Step ',num2str(step),': interpolation to spatial grid']));

definx = isfinite(Tomf); Tx=Tomf(definx);zxT=zTog(definx);
definx = isfinite(Lomf); Lx=Lomf(definx);zxL=zLog(definx);

% TEMPERATURE
zToi = z;
Toi = interp1(zxT,Tx,zToi,'cubic',NaN);
definT=isfinite(Toi);

zLoi = 0.5*(z(1:nz-1)+z(2:nz));
Loi = interp1(zxL,Lx,zLoi,'cubic',NaN);
definL=isfinite(Loi);
Lx=Loi(definL);zxL=zLoi(definL);
nn=length(Lx);Lmin=Lx(1);Lmax=Lx(nn);
Loi(zLoi <= min(zxL))=Lmin;
Loi(zLoi  > max(zxL))=Lmax;


if graph,
    % PLOT TEMPERATURES & GRADIENTS
    figure
    plot(To,zTo,'-b','LineWidth',1); hold on
    plot(Toi,zToi,'-r','LineWidth',2); hold on
    xlim([0 50]);
    ylim([0 2600]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['Outukumpu - Sept 5, 2008 (filtered,median \pm mad)']),'FontSize',fontsz,'FontWeight',fontwg);
    S1=strcat(['original data']);
    S2=strcat(['interpolated']);
    legend(S1,S2,'location', 'northeast')
    grid on
    file=strcat(['OKU2_temp.eps']);
    saveas(gcf,file,'epsc2')

    % PLOT CONDUCTIVITIY
    figure
    plot(Lo,zLo,'-b','LineWidth',1); hold on
    plot(Loi,zLoi,'-r','LineWidth',2); hold on

    ylim([0 2600]);
    xlim([0 10]);
    set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
    xlabel('\lambda (W m^{-1}K^{-1})','FontSize',fontsz,'FontWeight',fontwg)
    ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
    title(strcat(['Outukumpu - Sept 5, 2008 (filtered median \pm mad)']),'FontSize',fontsz,'FontWeight',fontwg)
    S1=strcat(['original data']);
    S2=strcat(['interpolated']);
    legend(S1,S2,'location', 'northeast')
    grid on
    file=strcat(['OKU2_Lamb.eps']);
    saveas(gcf,file,'epsc2')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP III: ESTIMATE Q0DepthInterval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if estq,
    step=step+1;
    disp(strcat([ 'Step ',num2str(step),': Estimate Q0 - only info']));
    for Q_EstDepth=zEst
        DepthInterval=zGrad>Q_EstDepth;
        Q = Qobs(DepthInterval);
        Q = Q(isfinite(Q));
        Q0m = median(Q);Q0t=mad(Q,1);
        disp(strcat(['Depth = ', num2str(Q_EstDepth) ' - 2500 m :', ...
            ' Qmed = ',num2str(Q0m*1000), ...
            ' Qmad = ',num2str(Q0t*1000), ...
            ' mW/m**2' ]));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP IV: DEFINE STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=step+1;
disp(strcat([ 'Step ',num2str(step),': setup model & data structures']));

% STORE INTO STRUCTURE <MODEL>
nzToi=length(zToi);
ip(1:nzToi-1)=[1:nzToi-1];
nip=length(ip);
nones=ones(nip,1)';
ip=ip';
k=Loi';
z=zToi;
nz=nzToi;
dz=diff(zToi);
gt=GST0;
qb=Q0;
p=por*nones';
r=2800*nones';
c=740*nones';
rc=r.*c;
h=0.*nones';
kA=0.0013*nones';
kB=0.0029*nones';

model.ip=ip;
model.k=k';
model.z=z;
model.nz=nz;
model.dz=dz;
model.gt=gt;
model.qb=qb;
model.p=p;
model.r=r;
model.c=c';
model.rc=rc;
model.h=0.*nones';
model.kA=0.0013*nones';
model.kB=0.0029*nones';

% STORE INTO STRUCTURE DATA
nz = length(zToi);
id=[1:nz];
id=find(definT);
nd=length(id);
Tobs=Toi;
zd=zToi;
cov=err^2*ones(1,nd);

data.Tobs=Tobs;
data.id=id;
data.nd=nd;
data.z=zd;
data.err=err;
data.cov=cov;

% SAVE
filename=strcat(['OKU_Qb',num2str(1000*Q0),'_final']);
save(filename);
disp([' ']);
disp(strcat([' >>>>> saved to ', filename]));
disp([' ']);

clear model data

% STRUCTURE MODEL
model= put_sitepar(k,kA,kB,h,p,c,r,rc,z,ip,qb,gt,site);
% STRUCTURE DATA
data= put_sitedat(Tobs,id,zd,cov,err,site);

file=strcat([tit '_site.mat']);
save(file,'data','model');
disp([' >>>>> site data saved to:' file]);
end
