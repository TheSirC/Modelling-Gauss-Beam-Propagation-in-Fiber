function beam=iso11146(pdd)
%ISO11146 Beam propagation ratios according to ISO 11146-1
%   BEAM=ISO11146(PDD) returns a strcut BEAM containing the propagation 
%   parameters according to the international standard ISO11146-1 "Lasers 
%   and laser-related equipment - test methods for laser beam widths, 
%   divergence angles and beam propagation ratios - part 1: Stigmatic and 
%   simple astigmatic beams" (2005), with the following (optional) fields
%
%   name       (Identifier, string)
%   filename   (Filename, string) 
%   wavelength (Laser wavelength [m], double)
%   z          (z-coordinates along optical axis, vector of doubles)
%   mx         (1st moment at all planes in x direction, vector of doubles)
%   my         (1st moment at all planes in y direction, vector of doubles)
%   vx         (2nd moment (variance) at all planes in x direction, vector of doubles)
%   vy         (2nd moment (variance) at all planes in y direction, vector of doubles)
%   dx         (Width of PDD at all planes in x direction, vector of doubles)
%   dy         (Length of PDD at all planes in x direction, vector of doubles)
%   dr         (Diameter of PDD at all planes in x direction, vector of doubles)
%   eta        (Elipticity dx/dy of PDD at all planes in x direction, vector of doubles)
%
%   If PDD defines three planes or more the following fields are also
%   returned
%
%   Cx         (Coefficients of 2nd order polynomial fit of PDD widths, vector of 3 doubles)
%   z0x        (Location of focus/waist [m] along optical axis in XZ plane, double)
%   d0x        (Width of focus/waist [m] in XZ plane, double)
%   divx       (Full fair-field divergence angle [mrad] in XZ plane, double)
%   zRx        (Rayleigh length [m] in XZ plane, double)
%   M2x        (Times-limited-diffraction number (beam quality) in XZ plane, double)
%   Cy         (Coefficients of 2nd order polynomial fit of PDD widths, vector of 3 doubles)
%   z0y        (Location of focus/waist [m] along optical axis in YZ plane, double)
%   d0y        (Length of focus/waist [m] in YZ plane, double)
%   divy       (Full fair-field divergence angle [mrad] in YZ plane, double)
%   zRy        (Rayleigh length [m] in YZ plane, double)
%   M2y        (Times-limited-diffraction number (beam quality) in YZ plane, double)
%   Cr         (Coefficients of 2nd order polynomial fit of PDD widths, vector of 3 doubles)
%   z0r        (Location of focus/waist [m] along optical axis, double)
%   d0r        (Diameter of focus/waist [m], double)
%   divr       (Full fair-field divergence angle [mrad], double)
%   zRr        (Rayleigh length [m], double)
%   M2r        (Times-limited-diffraction number (beam quality), double)
%
%   EXAMPLE
%   pdd=mdfread('example.mdf')
%   beam=iso11146(pdd)
%   dispbeam(beam)
%   caustic(beam)
%
%   See also dispbeam, causfit

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Initialisation

n=max(size(pdd)); %Number of planes in pdd

if isfield(pdd,'name')
    beam.name=pdd(1).name;
end;

if isfield(pdd(1),'filename')
   beam.filename=pdd(1).filename;
end;

if isfield(pdd(1),'date')
   beam.date=pdd(1).date;
end;

if isfield(pdd(1),'time')
   beam.time=pdd(1).time;
end;

if isfield(pdd(1),'wavelength')
    beam.wavelength=pdd(1).wavelength;
end;

if isfield(pdd(1),'z')
    beam.z=[]; 
end;

beam.mx=[]; %First moment (centroid) of beam [m] in x-direction
beam.my=[]; %First moment (centroid) of beam [m] in y-direction
beam.vx=[]; %Second moment (variance) of beam [m^2] in x-direction
beam.vy=[]; %Second moment (variance) of beam [m^2] in y-direction
beam.vxy=[]; %Cross variance of beam [m] in xy-plane
beam.dx=[]; %Beam width [m] in xz-plane
beam.dy=[]; %Beam length [m] in yz-plane
beam.dr=[]; %Beam diameter [m]
beam.eta=[]; %Ellipcity

%Calculations according to ISO1146
for i=1:n
 [mx,my,vx,vy,vxy]=momentxy(pdd(i).x,pdd(i).y,pdd(i).ixy); %Variances
 phi=atan(2*vxy/(vx-vy))/2; %Azimuth angle of asitgmatic beam [rad]
 gamma=sign(vx-vy);
 c1=vx+vy; %helper variable
 c2=gamma*sqrt((vx-vy)^2+4*vxy^2); %helper variable

 if isfield(pdd(1),'z')
    beam.z=[beam.z pdd(i).z];
 end;
 
 dx=2*sqrt(2)*sqrt(c1+c2);
 dy=2*sqrt(2)*sqrt(c1-c2);
 dr=2*sqrt(2)*sqrt(c1);

 beam.mx=[beam.mx mx];
 beam.my=[beam.my my];
 beam.vx=[beam.vx vx];
 beam.vy=[beam.vy vy];
 beam.vxy=[beam.vxy vxy];

 beam.dx=[beam.dx dx];
 beam.dy=[beam.dy dy];
 beam.dr=[beam.dr dr];
 

 eta=dx/dy; %Ellipticity
  if eta>1
    eta=1/eta;
 end;
 beam.eta=[beam.eta eta];
end; %for

if n>=3
    beam=causfit(beam);
end,

