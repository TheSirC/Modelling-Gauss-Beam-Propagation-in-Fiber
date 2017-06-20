function pdd=t2surfsrc(mat,T,v)
%T2SURFSRC Surface power density profile based on desired temperature profile.
%   PDD=T2SURFSRC(MAT,T,V) returns a PDD struct containing a power density 
%   distribution, which, when absorbed at the surface (z=0) generates the 
%   2D temperature profile struct T in a semi-inifite material MAT, when the
%   PDD moves at velocity V [m/s] relative to the material. The x and y 
%   dimensions of PDD are defined by the x and y coordinates of PDD. 
%
%   EXAMPLE:
%   load texample T
%   plottemp(T) 
%   load C45
%   pdd=t2surfsrc(C45,T,50e-3)
%   figure
%   plotpdd(pdd)
%   T2=tsurfsrc(C45,pdd,50e-3)
%   plottemp(T2)
%
%   See also GAUSS, TEMPL, TEMMN, TOPHAT, RECTUNIF

%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 20-sep-2010


if  prod(size(T))>1 
    error(['Laser Toolbox: T2SURFSRC: the temperature profile should be ',...  
    'defined at one time instance only.']);
end;

if isfield(T,'t')
    if ~isinf(T.t)
    warning(['Laser Toolbox: T2SURFSRC: steady state conditions (t=Inf) ',...
        'are considered.']);
    end;
end;
    
if isfield(T,'z')
    if prod(size(T.z))>1
       warning(['Laser Toolbox: T2SURFSRC: the temperature profile is considered',...
           ' to be defined for the surface (z=0).']);
    end;
    if (prod(size(T.z))==1 && T.z~=0)
       warning(['Laser Toolbox: T2SURFSRC: the temperature profile is considered',...
           ' to be defined for the surface (z=0).']);
    end;
end;

T.z=0;
lenx=length(T.x);
leny=length(T.y);

Tmin=min(min(T.Txyz));
if Tmin>0
   warning(['Laser Toolbox: T2SURFSRC: minimum temperature',...
       '(T=' num2str(Tmin) ') is substracted from the input.']);
   T.Txyz=T.Txyz-Tmin;
end;

dx=(max(T.x)-min(T.x))/(lenx-1);
dy=(max(T.y)-min(T.y))/(leny-1);
T.x=T.x-dx/2;
T.y=T.y-dy/2;
dA=dx*dy;
sigma=abs(v)/2/mat.kappa;

%Calculation of approximation Wo
beta=atan(dx/dy);
w1=sqrt((dx^2+4*T.z^2*cos(beta)^2)/(sin(beta)^2));
w2=sqrt((dy^2+4*T.z^2*sin(beta)^2)/(cos(beta)^2));
W1=(w1+dx)/(w1-dx);
W2=(w2+dy)/(w2-dy); 
W3=dx*log(W1)+dy*log(W2)-4*T.z*(atan(w1/2*T.z)+atan(w2/2*T.z))+2*T.z*pi;
Wo=1/2/mat.K/dx/dy/pi*exp(-v*T.z/2/mat.kappa)*W3;

C1=1/(2*pi*mat.K);

for i=1:lenx
   for j=1:leny
      R=sqrt(T.x(i)^2+T.y(j)^2+T.z^2);
      if (C1*exp(-sigma*(T.x(i)+R)))==0
         invW=Inf;
      else
         invW=R/(C1*exp(-sigma*(T.x(i)+R)));
      end;
      if (invW*Wo)<1
         W(j,i)=Wo;
      else
         W(j,i)=1/invW;
      end;
   end;
end;

FW=fft2(W);
if ~isempty(find(FW==0))
   warning('Laser Toolbox: I2T: singularity in Fourier transform.');
%   gamma=(max(max(abs(FW))))*0.01;
%   gamma=0;
%   FW=FW+gamma;
end;

FT=fft2(T.Txyz); %Fourier transform of T
FIxy=(FT./FW); %Fourier transform of Ixy
Ixy=fftshift(real(ifft2(FIxy)))/dA/mat.A;
P=sum(sum(Ixy))*dA; %Total power [W]

pdd.name='T2SURFSRC';
pdd.date=date;
c=clock;
pdd.time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ]; 
pdd.x=T.x;
pdd.y=T.y;
pdd.z=T.z;
pdd.ixy=Ixy;
pdd.power=P;










