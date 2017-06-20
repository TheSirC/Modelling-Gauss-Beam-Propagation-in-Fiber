function [] = calcul_x()
%% 2D Simulation of the gaussian beam through an optical fiber
% Every parameter in this program is in SI units.
%% Parameters
% Physical parameters

global l; % Wavelength
global Ns; % Silica idex
global Na; % Air index (for later)
%%
% Laser parameters

global On; % Numerical aperture
global w0; % Waist
global M; % Square root of the beam-quality factor
global zr; % Rayleigh distance
%%
% Fiber parameters

global R; % Curvature radius
global zi; % Position of the interface
global zmax; % Radius of the fiber
%%
% Numerical parameters

global res; % Numerical resolution
global x_window_width;
x = linspace(-x_window_width*w0,x_window_width*w0,res); % Transverse coordinates
z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates

% Global inline functions 
global q;
global w;

%% Simulations
% Matrix representation of the picture

Px = [];

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
    
    %     % Computation of the real radius
    %     R = 1/(1/Q(1,idx)+1i*l*M^2/(pi*Ns*W(1,idx)^2));
    %     Pr = [Pr,R];
    
    % Computation of the intensity on the zTemp-plane
    U = @(x,zTemp) 1/Q(1,idx)^(M^2/2)... <- For readability
        .*(exp(-(1i*pi*x.^2)/(l*Q(1,idx)))...
        .*polyval(hermitePoly(M^2),sqrt(2)*(M.*x)/(W(1,idx))));
    I = abs(feval(U,x,zTemp).^2);
    Px = vertcat(Px,I);
end
%% Plotting
% Intensity

figure; imagesc(Px'); colormap(hot); colorbar;
ax = gca;
pIx = find(z >= zi); % Retreving the position of the interface in the matrix
pFx = find(diff(sign(diff(var(Px,0,2)))) == 2); % Retreving the position of focalisation in the matrix...
                                                % by calculating the first zero of the first derivative of the variance...
                                                % AFTER the interface

msgbox(sprintf('The position of focalisation is calculated to be at %gm after the interface.',... <- Displaying the position in a message box
    ((z(1,end)-z(1,1))*(pFx(2,1)-pIx(1,1))/res)));
ax.XTickLabelRotation = 45;
xlabel('z values'); ax.XTick = [0 pIx(1,1) pFx(2,1) res-1]; ax.XTickLabel = {'0','interface','focalisation','center'};
ylabel('x values'); ax.YTick = [1 res/2 res]; ax.YTickLabel = {sprintf('-%d\omega_0','0','%d\omega_0',x_window_width,x_window_width)};

end