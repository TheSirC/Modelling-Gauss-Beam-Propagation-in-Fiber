%% Simulation of the 3D contribution for producing animation
% Every parameter in this program is in SI units.

clear; close all;
answer = {'1000','1','1.45','0.4','3','62.5','-40','3200'};
defaultans = {'1000','1','1.45','0.4','3','62.5','-40','3200'};

%% Parameters
% Physical parameters

global l;
l = check_user_input(answer,1,defaultans) * 1e-9; % Wavelength
global Ns;
Ns = check_user_input(answer,3,defaultans); % Silica idex
global Na;
Na = check_user_input(answer,2,defaultans); % Air index (for later)
%% 
% Laser parameters

global On;
On = check_user_input(answer,4,defaultans); % Numerical aperture
global w0;
w0 = check_user_input(answer,5,defaultans) * 1e-6; % Waist
global M;
M = sqrt(pi*w0^2*On/l); % Square root of the beam-quality factor
global zr;
zr = w0/On; % Rayleigh distance
%% 
% Fiber parameters

global R;
R = check_user_input(answer,6,defaultans) * 1e-6; % Curvature radius
global zi;
mkdir('Résultats')
for zi = -50e-6:10e-6:-10e-6
    cd('Résultats')
    mkdir(sprintf('%0.5e µm',-zi))
    cd(sprintf('%0.5e µm',-zi))
    global zmax;
    zmax = zi + R; % Radius of the fiber
    %%
    % Numerical parameters
    
    global res;
    res = check_user_input(answer,8,defaultans); % Numerical resolution
    global window_width;
    window_width = 1;
    volume_window = 250;
    global z;
    z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates
    global pix2metersXY;
    pix2metersXY = 2*window_width*w0/res; % Conversion factor to pass from pixels measurements to SI units for X and Y axis
    global pix2metersZ;
    pix2metersZ = (z(1,end)-z(1,1))/res; % Conversion factor to pass from pixels measurements to SI units for Z axis
    
    %%
    % Inline functions
    
    global q;
    q = @(z) z - 1i*zr; % Inline function to compute the Complex beam parameter
    global w;
    w = @(z) w0*sqrt(1+(z/z(1,1)^2)); % Inline function to compute the Gaussian beam width
    
    %% Simulations
    % Computations for both axes
    [Px,pIx,pFx] = calcul_x;
    [Py,pIy,pFy] = calcul_y;
    
    % Non-linear parameter
    thr = 0.5e10; % Beam diameter
    nlthr = 2e10;
    
    % Computations for the ellipses
    figure;
    for ze = pFx(2,1)-volume_window:pFy(2,1)+volume_window
        %[dX,dY] = calcul_ellipse(Px,Py,ze,thr);
        [dXn,dYn] = calcul_ellipse(Px,Py,ze,nlthr);
        X = @(t) dX*cos(t);
        Y = @(t) dY*sin(t);
        
        Xnl = @(t) dXn*cos(t);
        Ynl = @(t) dYn*sin(t);
        
        Z = @(t) ze;
        %fplot3(X,Y,Z,'LineWidth',0.25); hold on
        fplot3(Xnl,Ynl,Z,'r-'); hold on
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
    xlabel('x values'); axes.XTick = [-dX 0 dX]; axes.XTickLabel = {'0',sprintf('%.2e',res*pix2metersXY/2),sprintf('%.2e',pix2metersXY*res)};
    ylabel('y values'); axes.YTick = [-dY 0 dY]; axes.YTickLabel = {'0',sprintf('%.2e',res*pix2metersXY/2),sprintf('%.2e',pix2metersXY*res)};
    zlabel('z values'); axes.ZTick = [ pFx(2,1)-50 (pFx(2,1)+pFy(2,1))/2 pFy(2,1)+50 ]; axes.ZTickLabel = {sprintf('%.2e',pix2metersZ*(pFx(2,1)-50)),'non-linear center',sprintf('%.2g',(z(1,end)-z(1,1))*(pFy(2,1)+50)/res)};
    axes.ZAxis.Exponent = -6;
    
    
    % Writing 3D view to specific file
    print('IsometricView','-dpng');
    cd ..
end
