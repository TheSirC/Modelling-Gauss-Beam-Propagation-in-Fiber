function caustic(beam)
%CAUSTIC Caustic plot.
%   CAUSTIC(BEAM) plots the beam dimensions (width, length, diameter) of 
%   BEAM struct along the optical axis z. If three or more planes are 
%   defined in BEAM also the corresponding 2n order polynomal fits are
%   plotted.
%
%   EXAMPLE
%   pdd=mdfread('example.mdf')
%   beam=iso11146(pdd)
%   caustic(beam)
%   dispbeam(beam)
%
%   See also causfit, plotpdd, plottemp

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

if nargin~=1
    if nargin~=1
       error('Laser Toolbox: CAUSTIC: too many or no input arguments.')
    end;
end;

n=length(beam.z);
N=10*n; %number of z locations along optical axis in fit
if n<2
    error('LaserToolBox:nothing to plot.');
end;

z=beam.z;
dx=beam.dx;
dy=beam.dy;
dr=beam.dr;

%Graph
set(gcf,'Name',beam(1).name);
plot(-dx/2,z,'bs',-dy/2,z,'gd',-dr/2,z,'ro',dx/2,z,'bs',dy/2,z,'gd',dr/2,z,'ro');
legend('xz plane', 'yz plane' , 'polar');
ylabel('z [m]');
xlabel('r, x, y [m]');
if isfield(beam(1),'device')
   title(beam(1).device);
end;

if ~isempty(beam.Cx) %If the fit exists
    zfit=linspace(min(z),max(z),N);
    Dxfit=sqrt(polyval(beam.Cx,zfit));
    Dyfit=sqrt(polyval(beam.Cy,zfit));    
    Drfit=sqrt(polyval(beam.Cr,zfit));    
    hold;
    plot(-Dxfit/2,zfit,'b--',Dxfit/2,zfit,'b--');
    plot(-Dyfit/2,zfit,'g-.',Dyfit/2,zfit,'g-.');
    plot(-Drfit/2,zfit,'r-',Drfit/2,zfit,'r-');
    hold;
end;

grid on;
axis tight;

