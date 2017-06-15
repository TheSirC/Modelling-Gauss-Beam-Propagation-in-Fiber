%% Simulating the position of the focalized beam given the displament of the waist
%
%% Parameters
% Physical parameters

l = 1000e-6; % Wavelength
Ns = 1.45; % Silica idex
Na = 1; % Air index (for later)
%% 
% Laser parameters

On = 4e-1; % Numerical aperture
w0 = 3e-6; % Waist
M = sqrt(pi*w0^2*On/l); % Square root of the beam-quality factor
zr = w0/On; % Rayleigh distance
%% 
% Fiber parameters

R = 62.5e-6; % Curvature radius
zi = 0:-1e-6:-R ; % Position of the interface

%% 
% Numerical parameters

res = 1e3; % Numerical resolution 
x = linspace(-2*w0,2*w0,res); % Transverse coordinates
%% Simulations
% Inline function to compute the Complex beam parameter

q = @(z) z + 1i*zr;
%% 
% Inline function to compute the Gaussian beam width

w = @(z) w0*sqrt(1+(z/z(1,1)^2));
%% 
% Matrix representation of the picture

Wf = zeros(size(zi));
for dw = 1:numel(zi)
    zmax = zi(1,dw) + R; % Radius of the fiber
    z = linspace(-2*zmax,2*zmax,res); % Propagation coordinates
    
    W = NaN(size(z)); % Matrix containing all the values for the radius of the beam
    
    Q = NaN(size(z)); % Matrix containing all the complex beam parameters value
    Q(1,1) = feval(q,z(1,1)); % Initialization of the first value
    passed = false; % Boolean checking that we only pass the interface once (Numerical simulation artefact)
    
    A = 1;
    P = [];
    for idx = 2:numel(z)
        zTemp = z(idx); % Retrieving value corresponding to the index (Parallelization artefact)
        
        % Interface case
        if zTemp >= zi(1,dw) && ~passed % <- Not passed
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
            
            %         Q2 = Masf*[Q(idx-1);1];
            %         Q(1,idx) = Q2(1,1)/Q2(2,1); % Retrieving the first line and normalizing
            %         passed = true;
            
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
        P = vertcat(P,I);
    end
    pI = find(z >= zi(dw));
    [~,pF] = max(var(P,0,2));
    Wf(1,dw) = ((z(1,end)-z(1,1))*(pF-pI(1,1))/res); 
end