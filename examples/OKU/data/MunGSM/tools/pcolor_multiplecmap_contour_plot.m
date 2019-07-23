%pcolor_multiplecmap_contour_plot:
%This function generates 3 pcolor plots with an overlayed contour outline and each plot
%has a different colormap using an image capturing routine.
%bool_contour is used to toggle the contour outline on and off: {1,0}
%shade specifies the shading type for the shade.m function
%Example:
%pcolor_multiplecmap_contour_plot(peaks(25), peaks(15), peaks(35), bone, jet, hot, 1, 'interp')
%written by Max Kaznady, University of Toronto Atmospheric Physics Group

function pcolor_multiplecmap_contour_plot(data1, data2, data3, cmap1, cmap2, cmap3, bool_contour, shade)

disp('GRAPHING');
figure
clf;

%WARNING: if you do not have OpenGL, then... this is where you might experience a crash
set(0,'DefaultFigureRenderer','opengl')

%Set units and directions
set(gcf, 'units', 'in');
set(gca, 'YDir', 'normal')

%Constants for plot positions
XLEFT=0.2;
XRIGHT=0.4;
YSPACE=0.07;
UPSHIFT=0.07;
HEIGHT=0.29;

%First Graph
subplot('position', [XLEFT, 2/3+UPSHIFT, 1-XRIGHT, HEIGHT-YSPACE]);
%create plot
disp('Time for first graph:');
tic

colormap(cmap1);
contour(data1);
cax=caxis;
handle=pcolor(data1);
caxis(cax);
hold on
if strcmp(shade, 'interp') | strcmp(shade, 'flat') | strcmp(shade, 'faceted')
    shading(shade);
else
    shading('flat');
end
if bool_contour    
    [cs,h]=contour(data1, 'linecolor', 'k');
    %caxis([min(cons1), max(cons1)]);
    hc = get(h,'Children');
    cdata = cell2mat(get(hc,'Cdata'));
    i = find(isnan(cdata));
    tcdata = cdata(i-1);
    set(hc(tcdata<0),'LineStyle','-.')
    set(hc(tcdata==0),'LineStyle','-','LineWidth',1.2);
    set(hc(tcdata>=0),'LineStyle','-')
end
[r c]=size(data1);
%bottom
line('xdata', [1:c], 'ydata', ones(1,r), 'linewidth', 2);
%top
line('xdata', [1:c], 'ydata', r*ones(1,r), 'linewidth', 2);
%left
line('xdata', ones(1, c), 'ydata', [1:r], 'linewidth', 2);
%right
line('xdata', c*ones(1, c), 'ydata', [1:r], 'linewidth', 2);
hold off
a = get(gca,{'Xlim','Ylim'});
hi = image(frame2im(getframe(gca)));
set(hi,{'Xdata','Ydata'},a)
set(gca,{'Xlim','Ylim'},a,'Ydir','normal', 'tickdir', 'out')
toc

%Second Graph
subplot('position', [XLEFT, 1/3+UPSHIFT, 1-XRIGHT, HEIGHT-YSPACE]);
disp('Time for second graph:')
tic

colormap(cmap2);
contour(data2);
cax=caxis;
handle=pcolor(data2);
caxis(cax);
hold on
if strcmp(shade, 'interp') | strcmp(shade, 'flat') | strcmp(shade, 'faceted')
    shading(shade);
else
    shading('flat');
end
if bool_contour    
    [cs,h]=contour(data2, 'linecolor', 'k');
    %caxis([min(cons1), max(cons1)]);
    hc = get(h,'Children');
    cdata = cell2mat(get(hc,'Cdata'));
    i = find(isnan(cdata));
    tcdata = cdata(i-1);
    set(hc(tcdata<0),'LineStyle','-.')
    set(hc(tcdata==0),'LineStyle','-','LineWidth',1.2);
    set(hc(tcdata>=0),'LineStyle','-')
end
[r c]=size(data2);
%bottom
line('xdata', [1:c], 'ydata', ones(1,r), 'linewidth', 2);
%top
line('xdata', [1:c], 'ydata', r*ones(1,r), 'linewidth', 2);
%left
line('xdata', ones(1, c), 'ydata', [1:r], 'linewidth', 2);
%right
line('xdata', c*ones(1, c), 'ydata', [1:r], 'linewidth', 2);
hold off
a = get(gca,{'Xlim','Ylim'});
hi = image(frame2im(getframe(gca)));
set(hi,{'Xdata','Ydata'},a)
set(gca,{'Xlim','Ylim'},a,'Ydir','normal', 'tickdir', 'out')
toc

%Third Graph
subplot('position', [XLEFT, UPSHIFT, 1-XRIGHT, HEIGHT-YSPACE]);
disp('Time for third graph:')
tic

colormap(cmap3);
contour(data3);
cax=caxis;
handle=pcolor(data3);
caxis(cax);
hold on
if strcmp(shade, 'interp') | strcmp(shade, 'flat') | strcmp(shade, 'faceted')
    shading(shade);
else
    shading('flat');
end
if bool_contour    
    [cs,h]=contour(data3, 'linecolor', 'k');
    %caxis([min(cons1), max(cons1)]);
    hc = get(h,'Children');
    cdata = cell2mat(get(hc,'Cdata'));
    i = find(isnan(cdata));
    tcdata = cdata(i-1);
    set(hc(tcdata<0),'LineStyle','-.')
    set(hc(tcdata==0),'LineStyle','-','LineWidth',1.2);
    set(hc(tcdata>=0),'LineStyle','-')
end
[r c]=size(data3);
%bottom
line('xdata', [1:c], 'ydata', ones(1,r), 'linewidth', 2);
%top
line('xdata', [1:c], 'ydata', r*ones(1,r), 'linewidth', 2);
%left
line('xdata', ones(1, c), 'ydata', [1:r], 'linewidth', 2);
%right
line('xdata', c*ones(1, c), 'ydata', [1:r], 'linewidth', 2);
hold off
a = get(gca,{'Xlim','Ylim'});
hi = image(frame2im(getframe(gca)));
set(hi,{'Xdata','Ydata'},a)
set(gca,{'Xlim','Ylim'},a,'Ydir','normal', 'tickdir', 'out')
toc
