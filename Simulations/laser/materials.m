%MATERIALS Sets and saves materials parameters.
%   MATERIALS defines structs with materials parameters and saves these 
%   to file. The structs contain the following fields
%   
%   name    (Indentifier, string)
%   A       (Absorptivity, ((vector of) double(s))
%   lambda  (Laser wavelength [m] for which A is defined, (vector of) double(s))
%   K       (Thermal conductivity [W/(m*K), double]
%   rho     (Density [kg/m^3], double)
%   Cp      (Thermal heat capcity [J/(kg*K)], double)
%   kappa   (Thermal diffusitivity [m^2/s], double)
%   T       (Temperature [K] for which the above parameters are defined, double)
%   Tm      (Melt temperature [K], double)
%   Tv      (Vaporization temperature [K], double)
%
%   EXAMPLE
%   materials
%   load C45
%   load Ti6Al4V

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente


%C45
C45.name='C45 High Grade Steel';
C45.A=0.7;
C45.lambda=10.6e-6;
C45.K=40; %Thermal conductivity [W/m/K]
C45.rho=8750; %Density [kg/m^3]
C45.Cp=690; %Specific heat capacitance [J/kg/K]
C45.kappa=C45.K/C45.rho/C45.Cp;
C45.Tm=1450+273;
C45.filename='C45.mat';
save C45.mat C45; 
disp([C45.name ' was saved to: C45.mat']);

%Ti6Al4V
Ti6Al4V.name='Ti6Al4V';
Ti6Al4V.A=1; %Value unknown therefore set to 1
Ti6Al4V.lambda=10.6e-6;
Ti6Al4V.K=6.8; %Thermal conductivity [W/m/K]
Ti6Al4V.rho=4428; %Density [kg/m^3]
Ti6Al4V.Cp=564; %Specific heat capacitance [J/kg/K]
Ti6Al4V.kappa=Ti6Al4V.K/Ti6Al4V.rho/Ti6Al4V.Cp;
Ti6Al4V.Tm=1650+273;
save Ti6Al4V.mat Ti6Al4V;
disp([Ti6Al4V.name ' was saved to Ti6Al4V.mat']);
