function y=laguerre(n,l,x)
%LAGUERRE Laguerre polynomial.
%   Y=LAGUERRE(N,L,X) returns the generalized Laguerre polynomial of order 
%   N and index L of X, where N and L are a non-negative integers and is X 
%   a real scalar.
%
%   EXAMPLE:
%   y=laguerre(1,2,3.3)
%
%   See also templ, hermite

%   G.R.B.E. Römer 13-dec-1994
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

if prod(size(n))>1
    error('Laser Toolbox: HERMITE: N should be a non-negative scalar integer.');
end;

if (n<0 || abs(round(n)-n)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: N should be a non-negative integer.');
end;

if prod(size(l))>1
    error('Laser Toolbox: HERMITE: L should be a non-negative scalar integer.');
end;

if (l<0 || abs(round(l)-l)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: L should be a non-negative integer.');
end;

if n==0
   y=1;
elseif n==1,
   y=-1*x+1+l;
elseif n>1
   y= (2*n+l-1-x)/n*laguerre(n-1,l,x) - (n+l-1)/n*laguerre(n-2,l,x);
end;

