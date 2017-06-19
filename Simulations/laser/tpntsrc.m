function T=tpntsrc(mat,P,v,x,y,z)
% TPNTSRC Temperature profile of point surface heat source.
%   T=TPNTSRC(MAT,P,V,X,Y) returns the 3D temperature profile at x, y and z
%   coordinates defined by vectors X, Y, and Z in a semi-inifite material 
%   MAT due to a point heat source of P [W] moving at velocity V [m/s]
%   relative to the material.
%
%   Example
%   x=linspace(-2e-3,2e-3,128)
%   y=x
%   z=linspace(0,1e-3,16)
%   load C45
%   T=tpntsrc(C45,10,10e-3,x,y,z)
%   plottemp(T)
%
%   See also tpntsrc, tsurfsrc


%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 20-sep-2010


epslion=1e-12;

%Make sure that x and y are row vectors
s=size(x);
if prod(s)>max(s)
   error('Laser Toolbox: TPNTSRC: x shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    x=x';
end;
s=size(y);
if prod(s)>max(s)
   error('Laser Toolbox: TPNTSRC: y shall be a row vector.');
end;
if s(1,1)>1 %If column vector
    y=y';
end;

C1=P/2/pi/mat.K;
C2=abs(v)/2/mat.kappa;

X=repmat(x,length(y),1);
Y=repmat(y',1,length(x));

T.name='TpntSrc';
T.date=date;
c=clock;
T.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ];
T.x=x;
T.y=y;
T.z=z;

%See Dowden, J.M., The mathemaics of thermal modeling; An introduction to
%the theory of laser material processing, Chapmann & Hall, 2003 and
%Römer, G.R.B.E., Modelling and control of laser surface treatment,
%PhD-thesis, University of Twente, Enschede, The Netherlands, 1999.

for m=1:length(z)
    R=sqrt(Y.^2+X.^2+z(m)^2);
    if min(min(R))==0
        R=R+epsilon;
    end
    T.Txyz(:,:,m)=C1*exp(-C2*(X+R))./R;
end;






