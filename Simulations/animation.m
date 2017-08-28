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
global q;
q = @(z) z + 1i*zr; % Inline function to compute the Complex beam parameter
global w;
w = @(z) w0*sqrt(1+(z/zr)^2); % Inline function to compute the Gaussian beam width

%% 
% Fiber parameters
results = {'Focalisation théorique','Focalisation calculée en x','Focalisation calculée en y','Résolution axiale','Distance inter-focus'}; 
global R;
R = check_user_input(answer,6,defaultans) * 1e-6; % Curvature radius
global zi;
mkdir('Résultats')
cd('Résultats')
for zi = -41e-6:1e-6:0
    close all;
    %% Simulations
    % Computations for both axes
    %% Parameters
    % Physical parameters
    
    global l; % Wavelength
    global Ns; % Silica idex
    global Na; % Air index (for later)
    
    %%
    % Laser parameters
    global M;
    global w0; % Waist
    global Pin;
    
    %%
    % Fiber parameters
    global zmax;    
    zmax = zi + R; % Radius of the fiber
    %%
    % Numerical parameters
    
    global res; % Numerical resolution
    global window_width;
    x = linspace(-window_width*w0,window_width*w0,res); % Transverse coordinates
    z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates
    global pix2metersZ;
    pix2metersZ = (z(1,end)-z(1,1))/res; % Conversion factor to pass from pixels measurements to SI units for Z axis
   
    % Global inline functions
    global q;
    global w;
    
    %% Simulations
    % Matrix representation of the picture
    Px = zeros(res-1,res);
    
    W = NaN(size(z)); % Matrix containing all the values for the radius of the beam
    
    Q = NaN(size(z)); % Matrix containing all the complex beam parameters value
    Q(1,1) = feval(q,z(1,1)); % Initialization of the first value
    passed = false; % Boolean checking that we only pass the interface once (Numerical simulation artefact)
    
    A = 1;
    
    for idx = 2:numel(z)
        zTemp = z(idx); % Retrieving value corresponding to the index
        
        % Interface case
        if zTemp >= zi && ~passed % <- Not passed the interface
            % Air/Silica curved interface ABCD matrix
            B = 0;
            C = (Na-Ns)/(R*Ns);
            D = Na/Ns;
            Mas = [A,B;
                C,D];
            
            % Computation of the simulated output complex beam parameter after the interface
            Q2 = Mas*[Q(idx-1);1];
            Q(1,idx) = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
            passed = true;
            
        else
            % Free space matrix (air or silica)
            B = zTemp - z(idx-1); % Position of the considered point (variable)
            C = 0;
            D = 1;
            Mp = [A,B;
                C,D];
            
            % Computation of the simulated output complex beam parameter
            Q2 = Mp*[Q(idx-1);1];
            Q(1,idx) = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
        end
        
        W(1,idx) = feval(w,zTemp); % Keeping the current waist to find the minimum at the end
        
        U = @(x,zTemp,idx) 1/Q(1,idx)^(1/2)... <- For readability
            .*(exp(-(1i*pi*x.^2)/(l*Q(1,idx)))); % Inline function to compute the amplitude field
        
        I = 1/2*abs(feval(U,x,z,idx).^2);
        Px(idx,:) = I;
    end
    %% Plotting
    % Intensity
    pIx = find(z >= zi); % Retreving the position of the interface in the matrix
    pFx = find(diff(sign(diff(var(Px,0,2)))) == -2); % Retreving the position of focalisation in the matrix...
    % by calculating the third zero of the first derivative of the variance...
    % AFTER the interface
    
    
    %% Simulations
    % Matrix representation of the picture
    Py = zeros(res-1,res);
    
    W = NaN(size(z)); % Matrix containing all the values for the radius of the beam
    
    Q = NaN(size(z)); % Matrix containing all the complex beam parameters value
    Q(1,1) = feval(q,z(1,1)); % Initialization of the first value
    passed = false; % Boolean checking that we only pass the interface once (Numerical simulation artefact)
    
    A = 1;
    
    for idx = 2:numel(z)
        zTemp = z(idx); % Retrieving value corresponding to the index
        
        % Interface case
        if zTemp >= zi && ~passed % <- Not passed the interface
            % Air/Silica flat interface ABCD matrix
            B = 0;
            C = 0;
            D = Na/Ns;
            Mas = [A,B;
                C,D];
            
            % Computation of the simulated output complex beam parameter after the interface
            Q2 = Mas*[Q(idx-1);1];
            Q(1,idx) = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
            passed = true;
            
        else
            % Free space matrix (air or silica)
            B = zTemp - z(idx-1); % Position of the considered point (variable)
            C = 0;
            D = 1;
            Mp = [A,B;
                C,D];
            
            % Computation of the simulated output complex beam parameter
            Q2 = Mp*[Q(idx-1);1];
            Q(1,idx) = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
        end
        
        W(1,idx) = feval(w,zTemp); % Keeping the current waist to find the minimum at the end
        
        U = @(y,zTemp,idx) 1/Q(1,idx)^(1/2)... <- For readability
            .*(exp(-(1i*pi*y.^2)/(l*Q(1,idx)))); % Inline function to compute the amplitude field
        
        I = 1/2*abs(feval(U,x,z,idx).^2);
        Py(idx,:) = I;
    end
    %% Plotting
    % Intensity
    pIy = find(z >= zi); % Retreving the position of the interface in the matrix
    pFy = find(diff(sign(diff(var(Py,0,2)))) == -2); % Retreving the position of focalisation in the matrix...
    % by calculating the first zero of the first derivative of the variance...
    % AFTER the interface
    
    % 
    results = [results;{zi,(pix2metersZ*(pFx(2,1)-pIx(1,1))),(pix2metersZ*(pFy(2,1)-pIy(1,1))),pix2metersZ*res,pix2metersZ*abs(pFy(2,1)-pFx(2,1))}]
end

cell2csv('results.csv',results);
