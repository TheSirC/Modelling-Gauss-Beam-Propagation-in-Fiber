function pdd=tophat(P,d,varargin)
%TOPHAT Top hat power density profile.
%   PDD=TOPHAT(P,d) returns a structure (struct) PDD with a Top hat 
%   power density distribution of power P [Watt] and diameter d [m].
%
%   PDD=TOPHAT(P,d,N) returns the struct PDD of which the Top hat 
%   power density PDD.ixy [W/m^2] is a N-by-N matrix, and the corresponding
%   x and y vectors have length N. If not specified N equals 128.
%
%   PDD=TOPHAT(P,d,X,Y) returns the struct PDD of the power density
%   distribution at locations defined by the vectors X and Y.
%
%   EXAMPLES:
%   pdd=tophat(1000,3e-3)
%   plotpdd(pdd)
%
%   pdd=tophat(1000,3e-3,256)
%   plotpdd(pdd)
%
%   x=linspace(-1.5e-3,1.5e-3,64)
%   y=linspace(-1.5e-3,1.5e-3,128)
%   pdd=tophat(1000,2e-3,x,y)
%   plotpdd(pdd)
%
%   See also TEMPL, TEMMN, GAUSS, RECTUNIF, T2SURFSRC

%   G.R.B.E. Römer 23-aug-208
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Initialize
Ndefault=128;
W=1.2; %Diameter "multiplier" to calculate x and y, if x and y where not defined

if (nargin<2)
     error('Laser Toolbox: TOPHAT: too few input arguments.');
end;
 
if (nargin>4)
     error('Laser Toolbox: TOP HAT: too many input arguments.');
end;

if prod(size(P))>1
    error('Laser Toolbox: TOP HAT: laser power P must be scalar.');
end;

if prod(size(d))>1
    error('Laser Toolbox: TOPHAT: beam diameter d must be scalar.');
end;

switch nargin
    case 2 %No N, nor x and y are specified
        N=Ndefault;
        pdd.x=linspace(-W*abs(d),W*abs(d),N);
        pdd.y=pdd.x;
    case 3 %N is specified, but not x, nor y       
        N=varargin{1};
        pdd.x=linspace(-W*abs(d),W*abs(d),N);
        pdd.y=pdd.x;
    case 4 %N is not specified, but x and y are
        x=varargin{1};
        y=varargin{2};
        if ( (length(x)<2) | (length(y)<2))
            error('Laser Toolbox: TOPHAT: x and y shall be vectors.');
        end;
        pdd.x=x;
        pdd.y=y;
end; %switch

pdd.name='TOPHAT';
pdd.date=date;
c=clock;
pdd.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];

pdd.ixy=zeros(length(pdd.y),length(pdd.x));
pdd.power=P;
Io=4*P/(pi*d^2);
for i=1:length(pdd.x)
	for j=1:length(pdd.y)
       R=sqrt(pdd.x(i)^2+pdd.y(j)^2);
       if  R<=d/2 
           pdd.ixy(j,i)=Io;
       end;
    end;
end;

