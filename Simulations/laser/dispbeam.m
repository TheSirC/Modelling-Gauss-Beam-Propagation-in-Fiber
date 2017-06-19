function cs=dispbeam(beam)
%DISPBEAM Display beam characteristics. 
%   DISPBEAM(BEAM) displays the characteristics of laser beam struct BEAM.
%   
%   S=DISPPDD(BEAM) returns a cell array of strings of the characteristics 
%   of the laser beam.
%
%   EXAMPLE
%   pdd=mdfread('example.mdf')
%   beam=iso11146(pdd)
%   dispbeam(beam)
%
%   See also disppdd

%   G.R.B.E. Römer 3-mar-2010
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

dig='%0.2E'; %format of data string (see SPRINTF for details)

if isfield(beam,'z')
    n=length(beam.z); %Number of planes
else
    n=1; 
end;

k=1;
if isfield(beam,'name')
    s{k}=['Name      : ' beam(1).name]; k=k+1;
end;
if (isfield(beam,'datestr') & isfield(beam,'timestr'))
   s{k}=['Date/time  : ' beam(1).datestr '/' beam(1).timestr]; k=k+1;
end;
if isfield(beam,'filename')
    s{k}=['File      : ' beam(1).filename]; k=k+1;
end;   
if isfield(beam,'fileid')
   s{k}=['File ID    : ' beam(1).fileid]; k=k+1;
end;
if isfield(beam,'comment')
   s{k}=['Comments   : ' beam(1).comment]; k=k+1;
end;
if isfield(beam,'wavelength')
   s{k}=['Wavelength: ' num2pstr(beam(1).wavelength) 'm']; k=k+1;
end;


if n>2 %So if there was a beam fit over 3 or more planes

    s{k}=['-----------------------------------------']; k=k+1;
    s{k}=['Beam propagation parameters (of beam fit)' ]; k=k+1;
    s{k}=['-----------------------------------------']; k=k+1;

    s{k}=['Polar coordinates' ]; k=k+1;
    if isfield(beam,'d0r')    
       s{k}=['  Focus diameter (d0r)   : ' num2pstr(beam(1).d0r) 'm']; k=k+1;
    end;
    if isfield(beam,'z0r')    
       s{k}=['  Focus location (z0r)   : ' num2pstr(beam(1).z0r) 'm']; k=k+1;
    end;
    if isfield(beam,'divr')    
       s{k}=['  Divergence angle (divr): ' num2pstr(beam(1).divr) 'rad']; k=k+1;
    end;
    if isfield(beam,'zRr')    
       s{k}=['  Rayleigh length (zRr)  : ' num2pstr(beam(1).zRr) 'm']; k=k+1;
    end;
    if isfield(beam,'M2r')    
       s{k}=['  Beam quality M2r       : ' num2str(beam(1).M2r) ]; k=k+1;
    end;
    if isfield(beam,'Cr')    
       s{k}=['  Fit coefficients       : ' num2str(beam(1).Cr) ]; k=k+1;
    end;

    s{k}=['XZ-plane' ]; k=k+1;
    if isfield(beam,'d0x')    
       s{k}=['  Focus diameter (d0x)   : ' num2pstr(beam(1).d0x) 'm']; k=k+1;
    end;
    if isfield(beam,'z0x')    
       s{k}=['  Focus location (z0x)   : ' num2pstr(beam(1).z0x) 'm']; k=k+1;
    end;
    if isfield(beam,'divx')    
       s{k}=['  Divergence angle (divx): ' num2pstr(beam(1).divx) 'rad']; k=k+1;
    end;
    if isfield(beam,'zRx')    
       s{k}=['  Rayleigh length (zRx)  : ' num2pstr(beam(1).zRx) 'm']; k=k+1;
    end;
    if isfield(beam,'M2x')    
       s{k}=['  Beam quality M2x       : ' num2str(beam(1).M2x) ]; k=k+1;
    end;
    if isfield(beam,'Cx')    
       s{k}=['  Fit coefficients       : ' num2str(beam(1).Cx) ]; k=k+1;
    end;


    s{k}=['YZ-plane' ]; k=k+1;
    if isfield(beam,'d0y')    
       s{k}=['  Focus diameter (d0y)   : ' num2pstr(beam(1).d0y) 'm']; k=k+1;
    end;
    if isfield(beam,'z0y')    
       s{k}=['  Focus location (z0y)   : ' num2pstr(beam(1).z0y) 'm']; k=k+1;
    end;
    if isfield(beam,'divy')    
       s{k}=['  Divergence angle (divy): ' num2pstr(beam(1).divy) 'rad']; k=k+1;
    end;
    if isfield(beam,'zRy')    
       s{k}=['  Rayleigh length (zRy)  : ' num2pstr(beam(1).zRy) 'm']; k=k+1;
    end;
    if isfield(beam,'M2y')    
       s{k}=['  Beam quality M2y       : ' num2str(beam(1).M2y) ]; k=k+1;
    end;
    if isfield(beam,'Cy')    
       s{k}=['  Fit coefficients       : ' num2str(beam(1).Cy) ]; k=k+1;
    end;

end; % if n>2

for j=1:n 
   if n>1
      s{k}=['-----------------------------------------']; k=k+1;
      s{k}=['Plane       : ' int2str(j)]; k=k+1;       
      s{k}=['-----------------------------------------']; k=k+1;      
   end;
   if isfield(beam,'z')  
      s{k}=['z            : ' num2pstr(beam.z(j)) 'm']; k=k+1;
   end; 
   if isfield(beam,'mx')  
      s{k}=['mx           : ' num2str(beam.mx(j)) 'm']; k=k+1;
   end;   
   if isfield(beam,'my')  
      s{k}=['my           : ' num2str(beam.my(j)) 'm']; k=k+1;
   end;      
   if isfield(beam,'vx')  
      s{k}=['vx           : ' num2str(beam.vx(j)) 'm']; k=k+1;
   end;         
   if isfield(beam,'vy')  
      s{k}=['vy           : ' num2str(beam.vy(j)) 'm']; k=k+1;
   end;     
   if isfield(beam,'vxy')  
      s{k}=['vxy          : ' num2str(beam.vxy(j)) 'm^2']; k=k+1;
   end;     
   if isfield(beam,'dx')  
      s{k}=['Width (dx)   : ' num2pstr(beam.dx(j)) 'm']; k=k+1;
   end;   
   if isfield(beam,'dy')  
      s{k}=['Length (dy)  : ' num2pstr(beam.dy(j)) 'm']; k=k+1;
   end;  
   if isfield(beam,'dr')  
      s{k}=['Diameter (dr): ' num2pstr(beam.dy(j)) 'm']; k=k+1;
   end;     
end; %for





if nargout<1 %If no output parameter is specified write to stdout
   for k=1:max(size(s))
       disp(s{k});
   end;
else
    cs=s;
end; %if



