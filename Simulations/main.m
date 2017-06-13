clear all; close all;

%% Parameters
res = 200; % Numerical resolution 
l = 1000e-6; % Wavelength
On = 4e-1; % Numerical aperture
Ns = 1.45; % Silica idex
Na = 1; % Air index (for later)
R = 1e18; % Curvature radius
w0 = 1e-3; % Waist
M = sqrt(pi*w0^2*On/l); % Square root of the beam-quality factor
zr = w0/On; % Rayleigh distance
z0 = 1e2; % Position of the interface
zc = 50e-6; % Position of the considered point
zmax = 60e-6; % Radius of the fiber

%% Simulations
% Complex beam parameter
q = @(z) z - z0 + 1i*zr; % Inline function to compute the Complex beam parameter
% Computation of the second-order intensity moment
W = @(z) w0*sqrt(1+(z/z0)^2);

% Air/Silica interface ABCD matrix
A = 1;
B = 0;
C = (Na-Ns)/(R*Ns);
D = Na/Ns;
Mas = [A,B;
       C,D];

P = cell(res);
x = linspace(-3*w0,3*w0,res); % Transverse coordonates
for z = linspace(z0,zmax,res)
    
        % Propagation in the silica ABCD matrix
        B = z; % Position of the considered point (variable)
        C = 0;
        D = 1;
        Mp = [A,B;
            C,D];
        
        % Computation of the simulated output complex beam parameter at the point
        % considered
        Q2 = Mp*Mas*[q(z);1];
        q2 = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
        
        % Computation of the real radius
        % R = 1/(1/q2+1i*l*M^2/(pi*Ns*W^2));
        
        % Computation of the intensity on the z-plane
        U = @(x,z) 1/q2^(M^2/2)... <- For readability
                    .*(exp(-(1i*pi*x.^2)/(l*q2))...
                    .*polyval(hermitePoly(M^2),sqrt(2)*(M.*x)/(W(z))));
        I = abs(U(x,z).^2);
        P = [P;I];
end

%% Plotting
imagesc(P);