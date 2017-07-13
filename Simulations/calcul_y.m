function [Py,pIy,pFy] = calcul_y()
%% 2D Simulation of the gaussian beam through an optical fiber
% Every parameter in this program is in SI units.
%% Parameters
% Physical parameters

global l; % Wavelength
global Ns; % Silica idex
global Na; % Air index (for later)

%%
% Laser parameters

global w0; % Waist

%%
% Fiber parameters

global zi; % Position of the interface

%%
% Numerical parameters

global res; % Numerical resolution
global window_width;
y = linspace(-window_width*w0,window_width*w0,res); % Transverse coordinates
global z;
global pix2metersXY;
global pix2metersZ;

% Global inline functions 
global q;
global w;

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
    
    U = @(x,idx) sqrt(2/pi)*... % <- For readability purpose
        ((2*pi/l)*w0)/(2*pi*Q(1,idx))*...
        exp(1i*(2*pi/l)*(y.^2)/(2*Q(1,idx))); % Inline function to compute the amplitude field
    
    I = 1/2*abs(feval(U,y,idx).^2);
    Py(idx,:) = I;
end
%% Plotting
% Intensity
subplot(2,1,2); imagesc(Py'); colormap(hot); colorbar; hold off;
ay = gca;
pIy = find(z >= zi); % Retreving the position of the interface in the matrix
pFy = find(diff(sign(diff(var(Py,0,2)))) == -2); % Retreving the position of focalisation in the matrix...
                                                % by calculating the first zero of the first derivative of the variance...
                                                % AFTER the interface

msgbox(sprintf('The position of focalisation for the flat direction is calculated to be at %gm after the interface.',... <- Displaying the position in a message box
                (pix2metersZ*(pFy(2,1)-pIy(1,1)))),'Success','Help'); 
ay.XTickLabelRotation = 45;
xlabel('z values'); ay.XTick = [0 pIy(1,1) pFy(2,1) res-1]; ay.XTickLabel = {'0','interface','focalisation','center'};
ylabel('y values'); ylabel('x values'); ay.YTick = [1 res/2 res]; ay.YTickLabel = {'0',sprintf('%.2e',res*pix2metersXY/2),sprintf('%.2e',pix2metersXY*res)};

% Writing images to specific file
print('PlanYZ','-dpng');


end