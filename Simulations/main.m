%% Simulation of the 3D contribution
% Every parameter in this program is in SI units.

clear; close all;
%% User interaction
prompt = {'Enter the wavelength \lambda of the laser (in nm) :',...
          'Enter the outside index (air for example) n_1 :',...
          'Enter the inside index (silica for example) n_2 :',...
          'Enter the numerical aperture of the objective used to focus the beam \theta :',...
          'Enter the waist w_0 (in \mu m) :',...
          'Enter the curvature of the fiber (its radius) R_f (in \mu m) :',...
          'Enter the focus position relative to the interface z_{focus} (in \mu m) :',...
          'Enter the number of points for the simulation :',...
          'Enter the energy per pulse of the laser (in nJ):', ...
          'Enter the duration of pulse of the laser (in fs):'};
dlg_title = 'Parameters for the simulations';
num_lines = 1;
defaultans = {'1023','1','1.45','0.4','3','62.5','-40','3200','200','150'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);

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
global Ep;
Ep = check_user_input(answer,9,defaultans) * 1e-9;
global tp;
tp = check_user_input(answer,10,defaultans) * 1e-15;
global Pin;
Pin = Ep/tp * sqrt(pi/2);
%% 
% Fiber parameters

global R;
R = check_user_input(answer,6,defaultans) * 1e-6; % Curvature radius
global zi;
zi = check_user_input(answer,7,defaultans) * 1e-6;  % Position of the interface
global zmax;    
zmax = zi + R; % Radius of the fiber
%% 
% Numerical parameters

global res;
res = check_user_input(answer,8,defaultans); % Numerical resolution
global window_width;
window_width = 5;
window_width_3D = 250;
step_points_3D = 100;
global z;
z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates
global pix2metersXY;
pix2metersXY = 2*window_width*w0/res; % Conversion factor to pass from pixels measurements to SI units for X and Y axis
global pix2metersZ;
pix2metersZ = (z(1,end)-z(1,1))/res; % Conversion factor to pass from pixels measurements to SI units for Z axis

%% 
% Inline functions

global q;
q = @(z) z + 1i*zr; % Inline function to compute the Complex beam parameter
global w;
w = @(z) w0*sqrt(1+(z/zr)^2); % Inline function to compute the Gaussian beam width

%% Simulations
% Computations for both axes
[Px,pIx,pFx] = calcul_x;
[Py,pIy,pFy] = calcul_y;

% Cmputations for the volume
calcul_3D(Px,pFx,Py,pFy,1.82e6,window_width_3D,step_points_3D);
