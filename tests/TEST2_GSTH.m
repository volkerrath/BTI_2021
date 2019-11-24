function [ierr]=SITE_Init(name)
% Site
% Prepare paeloclimate forcing

ierr = 0;

load('common.mat')

y2s=3600*24*365.25;s2y=1./y2s;

out= 1;
debug=0;%


plotit              = 1;

% CONSTANTS

site                = 'Test2';
gsth_form           = 'steps';
gsth_method         = 'linear';
gsth_file           = 'Test2_GSTH.dat';
gsth_smooth         = 0;
gst0                = 10.;
pom                 = -5.;

% SET PATHS

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);
addpath([strcat(['./local'])]);


if plotit
    
    set_graphpars
    tscal = 1.*s2y;
    tlimits([10 120000]);
    Tlimits[-10 10]);
    close all
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE SYN_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_GSTH_in.mat']);
if exist(F)
    disp([' ...',mfilename ' defaults overwritten from ', F])
    load(F); mstruct(gsth_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

disp(['   ']);
disp(strcat([ ' ... Generating GSTH for ', name]));


F=strcat([name,'_TimeGrid.mat']);
load(F)
nt = length(t);

Dgsth =importdata(gsth_file);    
tgsth   = Dgsth(:,1)*y2s;
Tgsth   = Dgsth(:,2);

if strcmp(gsth_form,'steps')
    

    Tgsth   = [Tgsth; Tgsth(end)];
    [Tgst]  = set_stpgst(t,tgsth,Tgsth,L,pom,0);
    
elseif strcmp(gsth_form,'points')
    
    [Tgst]  = set_pntgst(t,tgsth,Tgsth,method,debug);
    
    
else
        
        error([mfilename,': option <',gsth_form,'> not defined. EXIT']);
end

timeGSTH = t(:);
TempGsth = Tgst(:);
F=strcat([name,'_Init']);
disp([' ']);
disp([ 'results to ',F]);
save(F,'timeGSTH','TempGSTH');


if plotit
    figure
    plot(-t(:)*tscal,[Tgst(:)] ,'-b','LineWidth',3);hold on
    grid on;
    xlim(tlimits);
    ylim([Tlimits]);
    TXT=strrep(strcat([name]),'_',' ');
    textloc(TXT,'south','FontSize',0.5*fontsz,'FontWeight',fontwg);
    xlabel('Time BP/2000 (a)');
    ylabel('\Delta T (K)');
    
    set(gca,'xscale','lin','xdir','rev',...
        'FontSize',fontsz,'FontWeight',fontwg);
    S=strcat([name,'_GLinGSTH']);
    saveas(gcf,S,plotfmt)
    
    set(gca,'xscale','log','xdir','rev',...
        'xtick',[10 100 1000 10000 100000],...
        'FontSize',fontsz,'FontWeight',fontwg);
    S=strcat([name,'_GLogGSTH']);
    saveas(gcf,S,plotfmt);
end









