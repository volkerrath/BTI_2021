
clear all
close all
clc

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

% SET PATHS
pltpath='./';
datpath='./';
srcpath='../../../';
utlpath='../../../';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

% ONLY FOR PARRALLEL  EXECUTION
run_parallel=0;
parcors=   1;
plotit = 1;

save('common','srcpath','utlpath','datpath','pltpath','run_parallel','parcors'),


dfmt=1;ffmt='.zip';
%archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS
if plotit
    set_graphpars
    %plotfmt='epsc2';
    plotfmt='png';
end
%
site       = 'SYNB';
props       = 'syn';
prepstr       = '_QMesh';

name=[site prepstr];

Ntest=3;
a=25; b=35;
Qtest = -[a + (b-a).*rand(1,Ntest)]*2.3253*1e-3;

F=strcat([name,'_Delta50']);
load(F)


for Q=Qtest
    upper=find(Q>QB,1,'first');lower=upper-1;
    Qupper=QMesh(:,upper);
    Qlower=QMesh(:,lower);
    
    Qi=Qlower+(Q-QB(lower))/(QB(upper)-QB(lower)).*(Qupper-Qlower);
    if plotit
        figure
        plot([Qlower(:) Qi(:) Qupper(:)], zinit(:),'LineWidth',linwdt)
        set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
        grid on
        
        legend(['lower: ', num2str(-QB(lower)*1000,'%.4g')],...
            ['interp: ', num2str(-Q*1000,'%.4g')], ...
            ['upper: ', num2str(-QB(upper)*1000,'%.4g')],...
            'location', 'best')
        F=strcat([name,'_Delta50_',num2str(-Q*1000,'%.4g')]);
        saveas(gcf,F,plotfmt)
    end
end