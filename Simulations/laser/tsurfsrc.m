function T=tsurfsrc(mat,pdd,v,varargin)
%TSURFSRC Temperature profile of surface heat source.
%   T=TSURFSRC(MAT,PDD,V) returns steady state (time is infinity) 3D 
%   temperature profile T in a struct TEMP in a semi-inifite material MAT 
%   due to a surface heat source defined by a single plane in PDD moving at
%   velocity V [m/s] relative to the material. The x and y dimensions of T 
%   are defined by the x and y coordinates of PDD. The z coordinates (depth)
%   are 0, DH/2 and DH[m] by default where DH is the heat penetration depth.
%
%   T=TSURFSRC(MAT,PDD,V,TIME) returns an array of temperature structs T;
%   one struct per time instance defined by vector TIME. If not defined
%   TIME=Inf.
%
%   T=TSURFSRC(MAT,PDD,V,TIME,Z) returns an array of temperature structs T;
%   one struct per time instance defined by vector TIME, at z coordinates
%   [m] defined by vector Z.
%
%   T=TSURFSRC(MAT,PDD,V,TIME,Z,METHOD) returns an array of temperature 
%   structs T at time instances, at z coordinates [m] Z, uning the
%   calculation method METHOD. METHOD is a string equalling 'mulitint' 
%   (default), which slow(er) but accurate, or 'fft' which is fast but 
%   inacurrate.
%
%   EXAMPLES
%   pdd=temmn(0,1,10,[1e-3 2e-3],32)
%   load C45
%   T=tsurfsrc(C45,pdd,10e-3)
%   plottemp(T)
%
%   T=tsurfsrc(C45,pdd,10e-3,[1e-3 0.1 Inf])
%   plottemp(T)
%
%   T=tsurfsrc(C45,pdd,10e-3,[1e-3 0.1 Inf],[0 1e-3 2e-3])
%   plottemp(T)
%
%   pdd=temmn(0,1,10,[1e-3 2e-3],128)
%   T=tsurfsrc(C45,pdd,10e-3,[1e-3 0.1 Inf],[0 1e-3 2e-3],'fft')
%   plottemp(T)
%
%   See also tpntsrc, tlinesrc


%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 20-sep-2010


%T=tsrfsrc(mat,pdd,v,t,z,method)

if prod(size(pdd))>1
    error('Laser Toolbox: TSURFSRC: the power densty profile should contain 1 plane only.');
end;

%Initialisation
method='multiint'; %Calculation method 'multiint' or 'fft'
options.verbose=1;  %(1) messages on, (0) messages off
options.filter='No'; %Applies to 'fft' method

%Calculation of thermal penetration depth dh [m]
beam=iso11146(pdd);
dh=2*sqrt(mat.kappa*beam.dx/abs(v)); 

nargs=max(size(varargin));

switch nargs
    case 0 %So no additional parameters were specified
         t=Inf; %So quasi steady state model
         z=[0 dh/2 dh]; %Defines 3 planes IN the material
    case 1 %Only time (vector) is specified
         z=[0 dh/2 dh]; %Defines 3 planes IN the material
         t=varargin{1};
    case 2 %Time (vector) and depth (vector) are specified
         t=varargin{1};
         z=varargin{2};
         T.z=z; 
    case 3 %Time (vector), depth (vector) and method are specified
         t=varargin{1};
         z=varargin{2};
         T.z=z; 
         if isstr(varargin{3})
             method=varargin{3};
         else
             error('Laser Toolbox: TSURFSRC: method must be a string.');
         end;
        otherwise
        error('Laser Toolbox: TSURFSRC: incorrect number of input parameters.');
end; %switch

x=pdd.x;
y=pdd.y;

%Calculation of some constants
sigma=abs(v)/2/mat.kappa;
C1=1/2/pi/mat.K;

%Calculation dA (Delta area)
lenx=length(x);
leny=length(y);
lenz=length(z);
lent=length(t);
dx=(max(x)-min(x))/(lenx-1);%Delta x (assuming equidistant vector x)
dy=(max(y)-min(y))/(leny-1); %Delta y (assuming equidistant vector y)
dA=dx*dy; %Area dA
beta=atan(dx/dy);

if options.verbose
   disp('TSURFSRC: Please wait...');
   tic;
end;

      
switch lower(method)
    case 'multiint' %Multi integration
      xp=x; %x prime
      yp=y; %y prime
    for n=1:lent
      T(n).name='TSURFSRC';
      T(n).date=date;
      c=clock;
      T(n).time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ]; 
      T(n).x=x;
      T(n).y=y;
      T(n).z=z;
      T(n).t=t(n);
      for m=1:lenz
           if options.verbose
               disp(['t=' num2str(t(n)) ' [s] @ z=' num2str(z(m)) ' [m]']);
           end;
           for i=1:lenx
             for j=1:leny
              Ws=0;            
              for k=1:lenx
                  for l=1:leny
                      R=sqrt( (x(i)-xp(k))^2 + (y(j)-yp(l))^2 + z(m)^2);
                      invW=R/C1/exp(-sigma*(x(i)-xp(k)+R));
                      %Calculation of constants for approximation near singularity
                      w1=sqrt((dx^2+4*z(m)^2*cos(beta)^2)/(sin(beta)^2));
                      w2=sqrt((dy^2+4*z(m)^2*sin(beta)^2)/(cos(beta)^2));
                      W1=(w1+dx)/(w1-dx);
                      W2=(w2+dy)/(w2-dy); 
                      W3=dx*log(W1)+dy*log(W2)-4*z(m)*(atan(w1/2*z(m))+atan(w2/2*z(m)))+2*z(m)*pi;
                      W0=1/2/mat.K/dx/dy/pi*exp(-abs(v)*z(m)/2/mat.kappa)*W3;
                      if (invW*W0)<1 %So if (1-Wo*invW)>1
                          W=W0;
                      else
                          W=1/invW;
                      end;
                      Ws=Ws+W*pdd.ixy(l,k);
                  end; % for l
              end; %for k
              if isinf(t(n))
                 U=1;
              else
                 if t(n)==0
                   U=0;
                 else
                   R=sqrt(x(i)^2+y(j)^2+z(m)^2);
                   U=(1-erf((R-v*t(n))/2/sqrt(mat.kappa*t(n)))+exp(R*v/mat.kappa)*(1-erf((R+v*t(n))/2/sqrt(mat.kappa*t(n)))))/2;
              end;                 
              end;
              Tz(j,i)=Ws*U; %Temperature (rise) at z(m)
            end; %for j
           end; %for l
%           T(n).Txyz(:,:,m)=Tz*mat.A*dA;
             L2=(max(x)-min(x))*(max(y)-min(y));
             N2=lenx*leny; 
             T(n).Txyz(:,:,m)=Tz*mat.A/N2*L2;
             T(n).method='multiint';
      end; % for m
    end; %for n  

    

    case 'fft' %Fast Fourier Transform

    for n=1:lent
      T(n).name='TSURFSRC';
      T(n).date=date;
      c=clock;
      T(n).time=[int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6)) ]; 
      T(n).x=x;
      T(n).y=y;
      T(n).z=z;
      T(n).t=t(n);
      for m=1:lenz
        if options.verbose
           disp(['t=' num2str(t(n)) ' [s] @ z=' num2str(z(m)) ' [m]']);
        end;
        for i=1:lenx
         for j=1:leny
             R=sqrt(x(i)^2+y(j)^2+z(m)^2);
             invW=R/C1/exp(-sigma*(x(i)+R));
             %Calculation of constants for approximation near singularity
             w1=sqrt((dx^2+4*z(m)^2*cos(beta)^2)/(sin(beta)^2));
             w2=sqrt((dy^2+4*z(m)^2*sin(beta)^2)/(cos(beta)^2));
             W1=(w1+dx)/(w1-dx);
             W2=(w2+dy)/(w2-dy); 
             W3=dx*log(W1)+dy*log(W2)-4*z(m)*(atan(w1/2*z(m))+atan(w2/2*z(m)))+2*z(m)*pi;
             W0=1/2/mat.K/dx/dy/pi*exp(-abs(v)*z(m)/2/mat.kappa)*W3;
             if (invW*W0)<1 %So if (1-Wo*invW)>1
                 W(j,i)=W0; 
             else
                 W(j,i)=1/invW;
             end;
             
             %Time dependency
             
             if isinf(t(n))
               U=1;
             else
               if t(n)==0
                 U=0;
               else
                 R=sqrt(x(i)^2+y(j)^2+z(m)^2);
                 U=(1-erf((R-v*t(n))/2/sqrt(mat.kappa*t(n)))+exp(R*v/mat.kappa)*(1-erf((R+v*t(n))/2/sqrt(mat.kappa*t(n)))))/2;
               end;                 
             end;
             W(j,i)=W(j,i)*U;

         end; %for j
         end; %for i
        FW=fft2(W); %Fourier transform of W;
        if ~isempty(find(FW==0))
          if options.verbose
             warning('Laser Toolbox: TSURFSRC: singularity in Fourier transform.');
          end;
          FW=FW+max(max(real(FW)))/1e9;
        end;
        FI=fft2(pdd.ixy); %Fourier transform of Ixy [w/m^2];
        FT=(FW.*FI); %Fourier transform of T [K]
        T(n).Txyz(:,:,m)=fftshift(real(ifft2(FT)))*mat.A*dA;
        T(n).method='fft';
        end; %for m
      end; %for n 

    otherwise
      error(['Laser Toolbox: TSURFSRC: Method "' method '" unknown.']);
   end; %switch
       
    if options.verbose
       toc %Showing elapsed time
    end;

    

    

    
    
