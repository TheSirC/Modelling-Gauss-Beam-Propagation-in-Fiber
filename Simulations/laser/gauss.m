function pdd=gauss(P,d,varargin)
%GAUSS Gaussian power density profile.
%   PDD=GAUSS(P,d) returns a structure (struct) PDD with a Gaussian 
%   power density distribution of power P [Watt] and diameter d [m].
%
%   PDD=GAUSS(P,d,N) returns the struct PDD of which the Gaussian 
%   power density PDD.ixy [W/m^2] is a N-by-N matrix, and the corresponding
%   x and y vectors have length N. If not specified N equals 128.
%
%   PDD=GAUSS(P,d,X,Y) returns the struct PDD of a Gaussian power density
%   distribution at locations defined by the vectors X and Y.
%
%   EXAMPLES:
%   pdd=gauss(1000,3e-3)
%   plotpdd(pdd)
%
%   pdd=gauss(1000,3e-3,256)
%   plotpdd(pdd)
%
%   x=linspace(-1.5e-3,1.5e-3,64)
%   y=linspace(-1.5e-3,1.5e-3,128)
%   pdd=gauss(1000,3e-3,x,y)
%   plotpdd(pdd)
%
%   See also TEMPL, TEMMN, TOPHAT, RECTUNIF, T2SURFSRC

%   G.R.B.E. Römer 11-aug-2010
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Note: GAUSS is just a wrapper for TEMPL

nargs=max(size(varargin));
switch nargs
    case 0
        pdd=templ(0,0,P,d);
    case 1
        pdd=templ(0,0,P,d,varargin{1});
    case 2
        pdd=templ(0,0,P,d,varargin{1},varargin{2});
    otherwise
        error('Laser Toolbox: GAUSS: incorrect number of input arguments.');
end;

pdd.name='GAUSS'; %Overwrite name