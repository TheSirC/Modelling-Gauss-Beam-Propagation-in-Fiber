function beamo=causfit(beami)
%CAUSFIT Second order polynomial fit of beam dimensions.
%   BEAM=CAUSFIT(BEAM) returns a struct BEAM extended with beam propagation
%   ratios based on a second order fit of beam dimensions, according to the
%   international standard ISO11146-1 "Lasers and laser-related equipment-
%   test methods for laser beam widths, divergence angles and beam 
%   propagation ratios - part 1: Stigmatic and simple astigmatic beams" 
%   (2005).
%   
%   CAUSFIT is called by ISO11146(PDD) is PDD contains 3 planes or more.
%
%   See also iso11146

%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 03-March-2010

n=length(beami.z);

if n<3
    error('LaserToolbox: CAUSFIT:can not fit a caustic.');
end;
if n<5
    warning('LaserToolbox: CAUSFIT: No. of planes should be more than 5.');
end;

beamo=beami; %Copy the incomming beam struct to the outgoing beam struct

%Beam propagation parameters

%xz plane
beamo.Cx=polyfit(beamo.z,beamo.dx.^2,2);
a=beamo.Cx(3); b=beamo.Cx(2); c=beamo.Cx(1);
D=sqrt(4*a*c-b^2); %Helper variable
beamo.z0x=-b/2/c; %Location [m] of waist (or focus) along z-axis
beamo.d0x=1/2/sqrt(c)*D; %Focus diameter [m] in xz-plane
beamo.divx=sqrt(c); %Divergence
beamo.zRx=beamo.d0x/beamo.divx; %Rayleigh length
if isfield(beamo,'wavelength')
    M2x=pi/8/beamo.wavelength*D;
    if M2x>1
        beamo.M2x=M2x;
    end;
else
    warning('LaserToolbox: CAUSFIT: no wavelength defined, so M^2 can not be determined.');
end;

%yz plane
beamo.Cy=polyfit(beamo.z,beamo.dy.^2,2);
a=beamo.Cy(3); b=beamo.Cy(2); c=beamo.Cy(1);
D=sqrt(4*a*c-b^2);
beamo.z0y=-b/2/c; %location of waist (or focus) 
beamo.d0y=1/2/sqrt(c)*D; %diameter of waist (or focus)
beamo.divy=sqrt(c); %beam divergence 
beamo.zRy=beamo.d0y/beamo.divy; %Rayleigh length
if isfield(beamo,'wavelength')
    M2y=pi/8/beamo.wavelength*D;
    if M2y>1
        beamo.M2y=M2y;
    end;
else
    warning('LaserToolbox: CAUSFIT: no wavelength defined, so M^2 can not be determined.');
end;

%polar coordinates
beamo.Cr=polyfit(beamo.z,beamo.dr.^2,2);
a=beamo.Cr(3); b=beamo.Cr(2); c=beamo.Cr(1);
D=sqrt(4*a*c-b^2);
beamo.z0r=-b/2/c; %location of waist (or focus) 
beamo.d0r=1/2/sqrt(c)*D;  %diameter of waist (or focus)
beamo.divr=sqrt(c); %beam divergence 
beamo.zRr=beamo.d0r/beamo.divr; %Rayleigh length
if isfield(beamo,'wavelength')
    M2r=pi/8/beamo.wavelength*D;
    if M2r>1
        beamo.M2r=M2r;
    end;
else
    warning('LaserToolbox: CAUSFIT: no wavelength defined, so M^2 can not be determined.');
end;

