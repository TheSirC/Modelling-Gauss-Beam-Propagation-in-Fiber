function [mx,my,vx,vy,vxy]=momentxy(x,y,Ixy)
%MOMENTXY First and second order moments.
%   [MX,MY,VX,VY,VXY]=MOMENTSXY(X,Y,Ixy) returns the first order moments
%   MX [m] and MY [m] in x and y direction, of the power density profile 
%   Ixy, at the coordinates specified by vectors X and Y, as well as the
%   corresponding second order moments (variance) VX [m^2] and VY [m^2], 
%   and cross variance VXY.
%
%   MOMENTXY is called by ISO11146.
%
%   EXAMPLE
%   pdd=rectunif(100,[2e-3 1e-3])
%   [mx,my,vx,vy,vxy]-momentxy(pdd.x,pdd.y,pdd.ixy)
%
%   See also iso11146

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

p=sum(sum(Ixy)); %pseudo power

%Make sure that x and y are row vectors
s=size(x);
if prod(s)>max(s)
   error('Laser Toolbox: VARXYXY: x shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    x=x';
end;
s=size(y);
if prod(s)>max(s)
   error('Laser Toolbox: VARXYXY: y shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    y=y';
end;

%First moments
X=repmat(x,length(y),1);
mx=sum(sum(X.*Ixy))/p;
Y=repmat(y',1,length(x));
my=sum(sum(Y.*Ixy))/p;

%Second moments (variance)
vx=sum(sum((X-mx).^2.*Ixy))/p;
vy=sum(sum((Y-my).^2.*Ixy))/p;

%Cross variance
vxy=sum(sum((X-mx).*(Y-my).*Ixy))/p;

%DONE