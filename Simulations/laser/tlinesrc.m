function T=tlinesrc(mat,Q,v,x,y)
% TLINESRC Temperature profile of line heat source.
%   T=TLINESRC(MAT,Q,V,X,Y) returns the 2D temperature profile at x and y 
%   coordinates defined by vectors X and Y, in a semi-inifite material MAT 
%   due to a line heat source of Q [W/m] moving at velocity V [m/s]
%   relative to the material.
%
%   EXAMPLE
%   x=linspace(-2e-3,2e-3,128)
%   y=x
%   load C45
%   T=tlinesrc(C45,10,10e-3,x,y)
%   plottemp(T)
%
%   See also tpntsrc, tsurfsrc


%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 20-sep-2010

%Make sure that x and y are row vectors
s=size(x);
if prod(s)>max(s)
   error('Laser Toolbox: TLINESRC: x shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    x=x';
end;
s=size(y);
if prod(s)>max(s)
   error('Laser Toolbox: TLINESRC: y shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    y=y';
end;


T.name='TlineSrc';
T.date=date;
c=clock;
T.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];

C1=Q/2/pi/mat.K;
C2=v/2/mat.kappa;

%See Dowden, J.M., The mathemaics of thermal modeling; An introduction to
%the theory of laser material processing, Chapmann & Hall, 2003

X=repmat(x,length(y),1);
Y=repmat(y',1,length(x));

% X=repmat(x,length(y),1);
% Y=repmat(y,length(x),1);
R=sqrt(Y.^2+X.^2);
t=C1*exp(-C2*X).*besselk(0,C2*R); %Line source model

T.x=x;
T.y=y;
T.Txyz(:,:,1)=t;
