%% Simulation of the 3D contribution for producing animation
% Every parameter in this program is in SI units.

clear; close all;
answer = {'1023','1','1.45','0.4','3','62.5','-40','3200'};
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
% Numerical parameters

global res;
res = check_user_input(answer,8,defaultans); % Numerical resolution
global window_width;
window_width = 1;
window_width_3D = 250;
step_points_3D = 25;
global pix2metersXY;
pix2metersXY = 2*window_width*w0/res; % Conversion factor to pass from pixels measurements to SI units for X and Y axis
%% 
% Fiber parameters

global R;
R = check_user_input(answer,6,defaultans) * 1e-6; % Curvature radius
global zi;
mkdir('Résultats')
for zi = -41e-6:1e-6:-5e-6
    close all;
    cd('Résultats')
    mkdir(sprintf('%0.5e µm',-zi))
    cd(sprintf('%0.5e µm',-zi))
    global zmax;
    zmax = zi + R; % Radius of the fiber
     
    global z;
    z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates
  
    global pix2metersZ;
    pix2metersZ = (z(1,end)-z(1,1))/res; % Conversion factor to pass from pixels measurements to SI units for Z axis
    
    %%
    % Inline functions
    
    global q;
    q = @(z) z + 1i*zr; % Inline function to compute the Complex beam parameter
    global w;
    w = @(z) w0*sqrt(1+(z/z(1,1)^2)); % Inline function to compute the Gaussian beam width
    
    %% Simulations
    % Computations for both axes
    [Px,pIx,pFx] = calcul_x;
    [Py,pIy,pFy] = calcul_y;
    
    % Cmputations for the volume
    calcul_3D(Px,pFx,Py,pFy,1e4,window_width_3D,step_points_3D);
    cd ../..
end
