function y=hermite(n,x)
%HERMITE Hermite polynomial.
%   Y=HERMITE(N,X) returns the Nth degree Hermite polynomial of X, where 
%   N is a non-negative integer and X a real scalar.
%  
%   EXAMPLE:
%   y=hermite(1,1)
%
%   See also temmn, laguerre

%   G.R.B.E. Römer 13-dec-1994
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

if prod(size(n))>1
    error('Laser Toolbox: HERMITE: n should be a non-negative scalar integer.');
end;

if (n<0 || abs(round(n)-n)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: n should be a non-negative integer.');
end;

if n==0
   y=1;
elseif n==1,
   y=2*x;
elseif n>1
   y=2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x);
end;