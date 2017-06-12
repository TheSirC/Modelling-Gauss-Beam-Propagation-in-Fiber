%% Parameters
l = 1000e-6; % Wavelength
On = 1; % Numerical aperture
Ng = 1.45; % Glass idex
Na = 1; % Air index (for later)
R = 1e18; % Curvature radius
z0c = 1e-9; % Waist
zrc = 1; % Rayleigh lenght
zc = 4; % Position of the considered point

%% Simulations
% Complex beam parameter
q = @(z,z0,zr) z - z0 + 1i*zr; % Inline function to calculate the Complex beam parameter 

% Air/Silica interface ABCD matrix
A = 1;
B = 0;
C = (Na-Ng)/(R*Ng);
D = Na/Ng;
Mas = [A,B;
       C,D];
   
% Propagation in the silica ABCD matrix
B = z; % Position of the considered point (variable)
C = 0;
D = 1;
Mp = [A,B;
      C,D];
   
% Calculate the simulated output complex beam parameter
Q2 = Mp*Mas*[q(zc,z0c,zrc);1];
q2 = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing 
