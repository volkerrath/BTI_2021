% analysis and plots of MC results
%keep nsteps maxiter it t year2sec NAME
clear
close all
clc

load MC-Elk1-100000-24_results.mat

c=M(:,1);

index1=~isnan(c);

valid=c(index1);

index=~isnan(M);
M=M(index);
M=reshape(M,length(valid),nsteps+3);

avM=zeros(1,nsteps);
for j=1:nsteps

    jM=M(:,j);
    a=median(jM);

    avM(j)=a;

end
mod=M(:,1:nsteps-1);

trange=[-15:0.1:5];



% dens=zeros(length(trange),nsteps);
% 
% 
% for i=1:nsteps-1
% 
%     for j=2:length(trange)
%         for k=1:length(mod)
%             if mod(k,i)<=trange(j) & mod(k,i)>trange(j-1)
%                 dens(j,i)=dens(j,i)+1;
%             end
%         end
%     end
% 
% end
% 
% % plotting gray code colored models within defined grid
% pcolor(dens)
% view(2);
% colormap(flipud(gray(256)));
% shading interp
% %set(gca,'LineStyle','none');

% hold on
% plot(nsteps-1,avM,'LineWidth',0.5,'Color','r','LineStyle','-')


        qmc=M(:,nsteps+1);
        hist(qmc,[0.02:0.001:0.05])


for i=1:length(M)
    mod=M(i,1:nsteps);
    ty=t/year2sec;
    [X,Y]=stairs(-ty,mod(it));



%         plot(X,Y,'LineWidth',0.5,'Color','b','LineStyle','-');hold on;
%         set(gca,'XScale','log','XDir','reverse')
%         xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
%         ylim([-11 7]);
%
%
%         text(1000,-8.00,['Iterations:',num2str(maxiter)],'FontSize',14)
%
%         title([NAME ': final model'],'FontSize',14)
%         grid on;


 end




ty=t/year2sec;
[X,Y]=stairs(-ty,avM(it));
figure
plot(X,Y,'LineWidth',0.5,'Color','b','LineStyle','-');hold on;
set(gca,'XScale','log','XDir','reverse')
xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
ylim([-16 6]);


text(1000,-14.00,['Iterations:',num2str(maxiter)],'FontSize',14)

title([NAME ': final model'],'FontSize',14)
grid on;



