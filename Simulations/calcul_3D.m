function [  ] = calcul_3D(Px,pFx,Py,pFy,thr,window_width_3D,step_points_3D)
%calcul_3D Function handling the computation and the printing of the 3D
% model
%  Px and Py are the computed plane for the ortogonal plane of propagation
%  pFx and PFy are the respective positions of focalisation on each plane
%  thr is the threshold to consider that the ellipse is part of the volume
%  window_width_3D is the range before and after both points of
%  focalisation you have to consider in the computation
%  step_points_3D is the step for the drawing (the smaller the ellipse
%  they'll be)

global z;
global zi;
global res;
global pix2metersXY;
global pix2metersZ;

% Maximum semi-minor (respectively major) axis fro printing purpose
dXmax = 0;
dYmax = 0;

% Computations for the ellipses
figure;
for ze = pFx(2,1)-window_width_3D:step_points_3D:pFy(2,1)+window_width_3D
    
    [dX,dY] = calcul_ellipse(Px,Py,ze,thr);
    
    dXmax = max([dX,dXmax]);
    dYmax = max([dY,dYmax]);

    X = @(t) dX*cos(t);
    Y = @(t) dY*sin(t);
    
    Z = @(t) ze;
    
    fplot3(X,Y,Z,'LineWidth',0.25,'Color','r'); hold on
end
hold off;

% Title
title(sprintf('pos = %.2e m',-zi));

% Formatting the axes
axes = gca;
axes.XTickLabelRotation = -45; 
axes.YTickLabelRotation = 45;
axes.ZAxis.TickLabelFormat = '%.2e m';
set(gca,'zdir','reverse'); % Switching z-axis direction (top to bottom ascending)
xlabel('x values'); axes.XTick = [-dXmax 0 dXmax]; axes.XTickLabel = {'0',sprintf('%.2e',res*pix2metersXY/2),sprintf('%.2e',pix2metersXY*res)};
ylabel('y values'); axes.YTick = [-dYmax 0 dYmax]; axes.YTickLabel = {'0',sprintf('%.2e',res*pix2metersXY/2),sprintf('%.2e',pix2metersXY*res)};
zlabel('z values'); axes.ZTick = [ pFx(2,1)-50 (pFx(2,1)+pFy(2,1))/2 pFy(2,1)+50 ]; axes.ZTickLabel = {sprintf('%.2e',pix2metersZ*(pFx(2,1)-50)),'non-linear center',sprintf('%.2g',(z(1,end)-z(1,1))*(pFy(2,1)+50)/res)};
axes.ZAxis.Exponent = -6;

% Writing 3D view to specific file
print('IsometricView','-dpng');
end

