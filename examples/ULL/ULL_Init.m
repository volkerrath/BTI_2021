function [ierr]=ULL_Init(name)
% Site OTUKUMPU
% Prepare paeloclimate initial

ierr = 0;

load('common.mat')

y2s=3600*24*365.25;s2y=1./y2s;

out= 1;
debug=0;%


plotit              = 1;
plotmod             = 0;
zlimits             =[0 2000];

% CONSTANTS

site                = 'ULL';
props               = 'ull';
init_type           = 'equi';
% inital_type = 'periodic';
%
GST0                = 7.;
Qb                  = -0.051;
POM                 = -4.;


% NUMERICAL PARAMETERS
out            = 1;

F=strcat([name,'_NUMPar.mat']);
load(F);

prepstr= ['_' init_type];
%prepstr= [''];

GSTH_file='LGC_Ull.csv';
L=3;
initial_iter=50;

% SET PATHS

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['src/props/',props])]);
addpath([srcpath,filesep,strcat(['tools'])]);

% dfmt=1;ffmt='.zip';
% archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);

F=strcat([name,'_NUMPar.mat']);
load(F);
mstruct(numpar)

% load MODEL
disp(strcat([' load model for ' name ]));
F=strcat([datpath name '_SITEPar.mat']);
load(F);
mstruct(sitepar);

F=[name '_init_in'];
if exist(F)
    load(F);
    mstruct(init_in);
end


if plotit
    set_graphpars
    close all
end

disp(['   ']);
disp(strcat([ ' ... Generating initial values for ', name]));


nt = length(t);dt=diff(t);
nz = length(z);dz=diff(z);

switch lower(init_type)
    case {'e' 'equi' 'equilibrium'}
        Tin=GST0-POM;
        T0=heat1dns(k, kA, kB,h,p,Qb,Tin,dz,ip,maxitnl,tolnl,freeze,out);
        Tinit = T0;
        
        zinit =z;
        if plotit
            figure;
            plot(Tinit,zinit,'-r','LineWidth',3);hold on
            grid on
            ylim(zlimits);
            ylabel('z (m)','FontSize',fontsz);
            xlabel('T (C)','FontSize',fontsz);
            set(gca,'YDir','reverse');
            set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
            dtext= strcat(['Initial temperatures, equi, props = ',props]);
            textloc(dtext,'northeast','FontSize',fontsz-4,'FontWeight',fontwg);
            
            file=strcat([site prepstr, '_Equi']);
            saveas(gcf,file,plotfmt);
        end
        
        
    case {'p' 'prior' 'periodic'}
        % >>>>>>>>>>>>>>>>>>>> TIME STEPPING
         
        % SETUP FORCING
        PGSTH =load(GSTH_file);
        tGSTH=PGSTH(:,1)*y2s;
        TGSTH=PGSTH(:,2);TGSTH=[GST0; TGSTH];
        [Tgst] = set_stpgst(t,TGSTH,tGSTH,L,0);
        %
 
        Tit = [];
        for iter=1:initial_iter
            
            if iter==1
                Tin=Tgst(1)+POM;
                T0=heat1dns(k, kA, kB,h,p,Qb,Tin,dz,ip,maxitnl,tolnl,freeze,out);
                Tinit = T0;
                Tref=Tgst(length(Tgst));
                Tr=heat1dns(k, kA, kB,h,p,Qb,Tref,dz,ip,maxitnl,tolnl,freeze,out);
            end
            
            GST=Tgst;  
       
            [Tcalc,G,Q,K]=heat1dnt(k,kA,kB,h,rc,p,Qb,...
                dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,freeze,out);
                
            T0old=T0;T0=Tcalc(:,nt);
            Tit=[Tit T0(:)];
            disp(['iteration #',num2str(iter),': norm (T0-T0old)/T0)= ',...
                    num2str(norm((T0(:)-T0old(:))/T0(:),'inf'))]);
        end
     
        if plotit
            figure
            ty= t*s2y;
            tscal = 1.; % 1e-3;1 
            plot(-ty(:)*tscal,[Tgst(:);GST0] ,'-b','LineWidth',3);hold on
            grid on;
            xlim([20 150000]);
            ylim([-15 20]);
            TXT=strcat([name,'',props,' Uniform']);
            textloc(TXT,'south','FontSize',0.5*fontsz,'FontWeight',fontwg)
            xlabel('Time BP/2000 (a)');ylabel('\Delta T (K)');
            
            set(gca,'xscale','lin','xdir','rev',...
                'FontSize',fontsz,'FontWeight',fontwg);
            S=strcat([name,'_GLinGSTH']);
            saveas(gcf,S,plotfmt)
            
            set(gca,'xscale','log','xdir','rev',...
                'xtick',[10 100 1000 10000 100000],...
                'FontSize',fontsz,'FontWeight',fontwg);
            S=strcat([name,'_GLo2gGSTH']);
            saveas(gcf,S,plotfmt);
            
            
            
            figure;
            plot(Tinit,z,':k','LineWidth',2);hold on
            plot(Tit,z,'LineWidth',1);hold on
            plot(Tit(:,initial_iter),z,'-r','LineWidth',3);hold on
            ylim(zlimits);    grid on
            ylabel('z (m)','FontSize',fontsz);
            xlabel('T (C)','FontSize',fontsz);
            set(gca,'YDir','reverse');
            set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
            dtext= strcat(['Initial temperatures, props = ',props,', iterations: ',num2str(initial_iter,'%i')]);
            textloc(dtext,'northeast','FontSize',fontsz-4,'FontWeight',fontwg);
            file=strcat([name, '_Iterations']);
            saveas(gcf,file,plotfmt);
            
            
            figure;
            plot(Tinit,z,'-r','LineWidth',3);hold on
            grid on
            ylim(zlimits);
            ylabel('z (m)','FontSize',fontsz);
            xlabel('T (C)','FontSize',fontsz);
            set(gca,'YDir','reverse');
            set(gca,'FontSize',fontsz, 'FontWeight',fontwg);
            dtext= strcat(['Initial temperatures, final, props = ',props]);
            textloc(dtext,'northeast','FontSize',fontsz-4,'FontWeight',fontwg);
            
            file=strcat([name, '_Final']);
            saveas(gcf,file,plotfmt);
            
        end
    otherwise
        
        error([mfilename,': option <',init_type,'> not defined. EXIT']);
end


zinit=z;
Tinit=T0;
Tinit
F=strcat([name,'_Init']);
disp([' ']);
disp([ 'results to ',F]);
save(F,'zinit','Tinit','POM');





