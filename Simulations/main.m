%% Simulation of the 3D contribution
% Every parameter in this program is in SI units.

clear; close all;
%% User interaction
prompt = {'Enter the wavelength \lamda of the laser (in \mu m) :',...
          'Enter the outside index (air for example) n_1 :',...
          'Enter the inside index (silica for example) n_2 :',...
          'Enter the numerical aperture of the objective used to focus the beam \theta :',...
          'Enter the waist w_0 (in \mu \!m) :',...
          'Enter the curvature of the fiber (its radius) R_f (in \mu m) :',...
          'Enter the focus position relative to the interface z_{focus} (in \mu m) :',...
          'Enter the number of points for the simulation :'};
dlg_title = 'Parameters for the simulations';
num_lines = 1;
defaultans = {'1000','1','1.45','4e-1','3e-6','62.5e-6','-40e-6','1.5e3'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);

%% Parameters
% Physical parameters

global l;
l = check_user_input(answer{1}) * 10e-6; % Wavelength
global Ns;
Ns = check_user_input(answer{3}); % Silica idex
global Na;
Na = check_user_input(answer{2}); % Air index (for later)
%% 
% Laser parameters

global On;
On = check_user_input(answer{4}); % Numerical aperture
global w0;
w0 = check_user_input(answer{5}) * 10e-6; % Waist
global M;
M = sqrt(pi*w0^2*On/l); % Square root of the beam-quality factor
global zr;
zr = w0/On; % Rayleigh distance
%% 
% Fiber parameters

global R;
R = check_user_input(answer{6}) * 10e-6; % Curvature radius
global zi;
zi = check_user_input(answer{7}) * 10e-6;  % Position of the interface
global zmax;
zmax = zi + R; % Radius of the fiber
%% 
% Numerical parameters

global res;
res = check_user_input(answer{8}); % Numerical resolution
global x_window_width;
x_window_width = 64;
global y_window_width;
y_window_width = 64;
global z;
z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates


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

% Computations for the ellipses
for ze = pFx(2,1)-200:pFy(2,1)+200
    [dX,dY] = calcul_ellipse(Px,Py,ze);
    X = @(t) dX*cos(t);
    Y = @(t) dY*sin(t);
    Z = @(t) ze;
    fplot3(X,Y,Z); hold on
end
hold off;