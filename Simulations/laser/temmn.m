function pdd=temmn(m,n,P,d,varargin)
%TEMMN Gauss-Hermite mode.
%   PDD=TEMmn(m,n,P,d) returns a structure (struct) PDD with a laser 
%   power density distribution of a passive laser resonator in Gauss-
%   Hermite mode (TEMmn) with scalar mode numbers m and n, power P [Watt]
%   and where diameter d [m] is a 2-element vector with dimensions of the 
%   distribution in x and y direction.
%
%   PDD=TEMMN(m,n,P,d,N) returns the struct PDD of which the power density
%   PDD.ixy [W/m^2] is a N-by-N matrix, and the corresponding x- and y-
%   vectors have length N. If not specified N equals 128.
%
%   PDD=TEMMN(m,n,P,d,X,Y) returns the struct PDD of the power density
%   distribution at locations defined by the vectors X and Y.
%
%   EXAMPLES:
%   pdd=temmn(1,0,1000,[1 3]*1e-3)
%   plotpdd(pdd)
%
%   pdd=temmn(1,0,1000,[1 3]*1e-3,256)
%   plotpdd(pdd)
%
%   x=linspace(-1.5e-3,1.5e-3,64)
%   y=linspace(-3.5e-3,3.5e-3,128)
%   pdd=temmn(1,0,1000,[1 3]*1e-3,x,y)
%   plotpdd(pdd)
%
%   See also TEMPL, GAUSS, TOPHAT, RECTUNIF, T2SURFSRC, HERMITE

%   G.R.B.E. Römer 23-aug-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Initialize
Ndefault=128;
W=1.2; %Diameter "multiplier" to calculate x and y, if x and y where not defined

if (nargin<4)
     error('Laser Toolbox: TEMMN: too few input arguments.');
end;
 
if (nargin>6)
     error('Laser Toolbox: TEMMN: too many input arguments.');
end;

switch nargin
    case 4 %No N, nor x and y are specified
        N=Ndefault;
        pdd.x=linspace(-W*max(d),W*max(d),N);
        pdd.y=pdd.x;
    case 5 %N is specified, but not x, nor y       
        N=varargin{1};
        pdd.x=linspace(-W*max(d),W*max(d),N);
        pdd.y=pdd.x;
    case 6 %N is not specified, but x and y are
        x=varargin{1};
        y=varargin{2};
        if ( (length(x)<2) | (length(y)<2))
            error('Laser Toolbox: TEMMN: x and y shall be vectors.');
        end;
        pdd.x=x;
        pdd.y=y;
end; %switch

if prod(size(m))>1
    error('Laser Toolbox: TEMMN: mode number m must be scalar.');
end;

if prod(size(n))>1
    error('Laser Toolbox: TEMMN: mode number n must be scalar.');
end;

if prod(size(P))>1
    error('Laser Toolbox: TEMMN: laser power P must be scalar.');
end;

if prod(size(d))~=2
    error(['Laser Toolbox: TEMMN: beam diameter d must be a 2 element ',...
        'vector.']);
end;

pdd.name=['Gauss-Hermite TEM_{' int2str(m) ',' int2str(n) '}'];
pdd.date=date;
c=clock;
pdd.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];

dx=d(1,1);
dy=d(1,2);

Io=P*2^(3-m-n)*sqrt((2*m+1)*(2*n+1))/factorial(m)/factorial(n)/pi/dx/dy;

for i=1:length(pdd.x)
	for j=1:length(pdd.y)
         Hx=hermite(m,2*pdd.x(i)*sqrt(2)*sqrt(2*m+1)/dx);
	     Hy=hermite(n,2*pdd.y(j)*sqrt(2)*sqrt(2*n+1)/dy);
         pdd.ixy(j,i)=Io*(Hx*exp(-4*pdd.x(i)^2*(2*m+1)/dx^2))^2*(Hy*exp(-4*pdd.y(j)^2*(2*n+1)/dy^2))^2;
	end;
   end;
pdd.power=P;
pdd.ixy=real(pdd.ixy);


