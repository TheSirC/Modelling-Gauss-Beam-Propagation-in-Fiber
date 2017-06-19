function pdd=rectunif(P,d,varargin)
%RECTUNIF Rectangular uniform power density profile.
%   PDD=RECTUNIF(P,d) returns a structure (struct) PDD with a rectangular 
%   uniform power density distribution of power P [Watt] and where diameter
%   d [m] is a 2-element vector with dimensions of the distribution in x 
%   and y direction.
%
%   PDD=RECTUNIF(P,d,N) returns the struct PDD of which the power density
%   PDD.ixy [W/m^2] is a N-by-N matrix, and the corresponding x- and y-
%   vectors have length N. If not specified N equals 128.
%
%   PDD=RECTUNIF(P,d,X,Y) returns the struct PDD of the power density
%   distribution at locations defined by the vectors X and Y.
%
%   EXAMPLES:
%   pdd=rectunif(1000,[1 3]*1e-3)
%   plotpdd(pdd)
%
%   pdd=rectunif(1000,[1 3]*1e-3,256)
%   plotpdd(pdd)
%
%   x=linspace(-1.5e-3,1.5e-3,64)
%   y=linspace(-3.5e-3,3.5e-3,128)
%   pdd=rectunif(1000,[1 3]*1e-3,x,y)
%   plotpdd(pdd)
%
%   See also TEMPL, TEMMN, GAUSS, TOPHAT, T2SURFSRC

%   G.R.B.E. Römer 23-aug-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Initialize
Ndefault=128;
W=1.2; %Diameter "multiplier" to calculate x and y, if x and y where not defined

if (nargin<2)
     error('Laser Toolbox: RECTUNIF: too few input arguments.');
end;
 
if (nargin>4)
     error('Laser Toolbox: RECTUNIF: too many input arguments.');
end;

if length(d)~=2
       error('Laser Toolbox: RECTUNIF: the beam dimensions d shall be a 2 element vector.');
end;

switch nargin
    case 2 %No N, nor x and y are specified
        N=Ndefault;
        pdd.x=linspace(-W*max(d),W*max(d),N);
        pdd.y=pdd.x;
    case 3 %N is specified, but not x, nor y       
        N=varargin{1};
        pdd.x=linspace(-W*max(d),W*max(d),N);
        pdd.y=pdd.x;
    case 4 %N is not specified, but x and y are
        x=varargin{1};
        y=varargin{2};
        if ( (length(x)<2) | (length(y)<2))
            error('Laser Toolbox: RECTUNIF: x and y shall be vectors.');
        end;
        pdd.x=x;
        pdd.y=y;
end; %switch

pdd.name='RECTUNIF';
pdd.date=date;
c=clock;
pdd.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];

dx=d(1); dy=d(2);
pdd.ixy=zeros(length(pdd.y),length(pdd.x));
pdd.power=P;
Io=P/(dx*dy);
for i=1:length(pdd.x)
	for j=1:length(pdd.y)
       if ( (abs(pdd.x(i))<(dx/2)) && (abs(pdd.y(j))<(dy/2)) )
           pdd.ixy(j,i)=Io;
       end;
    end;
end;

