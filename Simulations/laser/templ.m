function pdd=templ(p,l,P,d,varargin)
%TEMPL Gauss-Laguerre mode.
%   PDD=TEMPL(p,l,P,d) returns a structure (struct) PDD with a laser 
%   power density distribution of a passive laser resonator in Gauss-
%   Laguerre mode (TEMpl) with scalar mode numbers p and l, power P [Watt]
%   and diameter d [m].
%
%   PDD=TEMPL(p,l,P,d,N) returns the struct PDD of which the power density
%   PDD.ixy [W/m^2] is a N-by-N matrix, and the corresponding x- and y-
%   vectors have length N. If not specified N equals 128.
%
%   PDD=TEMPL(p,l,P,d,X,Y) returns the struct PDD of the power density
%   distribution at locations defined by the vectors X and Y.
%
%   EXAMPLES:
%   pdd=templ(1,0,1000,3e-3)
%   plotpdd(pdd)
%
%   pdd=templ(1,0,1000,3e-3,256)
%   plotpdd(pdd)
%
%   x=linspace(-1.5e-3,1.5e-3,64)
%   y=linspace(-1.5e-3,1.5e-3,128)
%   pdd=templ(1,0,1000,3e-3,x,y)
%   plotpdd(pdd)
%
%   See also TEMMN, GAUSS, TOPHAT, RECTUNIF, T2SURFSRC, LAGUERRE

%   G.R.B.E. Römer 23-aug-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Initialize
Ndefault=128;
W=1.2; %Diameter "multiplier" to calculate x and y, 
       %if x and y were not defined

if (nargin<4)
     error('Laser Toolbox: TEMPL: too few input arguments.');
end;
 
if (nargin>6)
     error('Laser Toolbox: TEMPL: too many input arguments.');
end;

if prod(size(p))>1
    error('Laser Toolbox: TEMPL: mode number p must be scalar.');
end;

if prod(size(l))>1
    error('Laser Toolbox: TEMPL: mode number l must be scalar.');
end;

if prod(size(P))>1
    error('Laser Toolbox: TEMPL: laser power P must be scalar.');
end;

if prod(size(d))>1
    error('Laser Toolbox: TEMPL: beam diameter d must be scalar.');
end;

switch nargin
    case 4 %No N, nor x and y are specified
        N=Ndefault;
        pdd.x=linspace(-W*d,W*d,N);
        pdd.y=pdd.x;
    case 5 %N is specified, but not x, nor y       
        N=varargin{1};
        pdd.x=linspace(-W*d,W*d,N);
        pdd.y=pdd.x;
    case 6 %N is not specified, but x and y are
        x=varargin{1}; %x vector
        y=varargin{2}; %y vector
        if ( (length(x)<2) | (length(y)<2))
            error('Laser Toolbox: TEMPL: x and y shall be vectors.');
        end;
        pdd.x=x;
        pdd.y=y;
end; %switch

pdd.name=['Gauss-Laguerre TEM_{' int2str(p) ',' int2str(l) '}'];
pdd.datestr=date;
c=clock;
pdd.timestr=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];

M2=2*p+l+1;
  if l==0
     Io=8*(2*p+l+1)*P/pi/d^2;
  else
     Io=16*(2*p+l+1)*factorial(p)*P/pi/d^2/factorial(p+l);
  end;
  for i=1:length(pdd.x)
     for j=1:length(pdd.y)
        r=sqrt(pdd.x(i)^2+pdd.y(j)^2);
        if r==0
           phi=sign(pdd.x(i))*pi/2;
        else
   	      phi=acos(pdd.x(i)/r);
        end;
        rho=8*r^2*M2/d^2;
        L=laguerre(p,l,rho);
        pdd.ixy(j,i)=Io*rho^l*L^2*(cos(l*phi))^2*exp(-rho);
      end; 
   end;
pdd.power=P;
pdd.ixy=real(pdd.ixy);


