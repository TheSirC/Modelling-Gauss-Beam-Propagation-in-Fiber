function O=overlap(d,f,v)
%OVERLAP Overlap plot of pulses
%   OVERLAP(D,F,V) plots 10 circular laser pulses of diameter D [m], at 
%   pulse frequency F [Hz], at velocity V [m/s] of the laser beam relative 
%   to the substrate.
%
%   O=OVERLAP(D,F,V) returns the overlap O [0-1] only.
%
%   EXAMPLES:
%   overlap(20e-6,200e3,0.8)
%   O=overlap(20e-6,200e3,0.8)
%    

%   G.R.B.E. Römer 20-sep-2010
%   Copyright 2010 G.R.B.E. Römer, University of Twente

N=10; %Number of spots drawn in graph when no output arguments are specified
n=25; %Number of points on half a circle (spot)
C='red'; % Color of circle (spot)

%Calculating overlap
OL=1-v/f/d;
if OL<0
    OL=0;
end

%Plotting

if nargout~=1
figure; hold;
   theta=linspace(0,2*pi,n);
   x=0.5*d*cos(theta);
   y=0.5*d*sin(theta);
   delta=v/f;
   for i=0:N-1
       fill(x,y+i*delta,C);
   end;
hold;

title(['Overlap=' num2str(OL*100) '%, d=' num2pstr(d) 'm, f=' num2pstr(f) 'Hz, ' ...
    'v=' num2pstr(v) 'm/s \Deltax=' num2pstr(v/f) 'm.']);
axis equal;
grid on;
xlabel('y [m]');
ylabel('x [m]');
else
    O=OL;
end; %nargout






