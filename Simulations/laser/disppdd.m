function cs=disppdd(pdd)
%DISPPDD Display power density distribution. 
%   DISPPDD(PDD) displays the characteristics of each plane of a power
%   density profile struct PDD.
%   
%   S=DISPPDD(PDD) returns a cell array of strings of the characteristics 
%   of each plane.
%
%   EXAMPLE
%   pdd=gauss(100,1e-3)
%   disppdd(pdd)
%
%   See also dispbeam

%   G.R.B.E. Römer 3-mar-2010
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

n=max(size(pdd));
dig='%0.2E'; %format of data string (see SPRINTF for details)

k=1;
if isfield(pdd,'name')
  s{k}=['Name     : ' pdd(1).name]; k=k+1;
end;
if (isfield(pdd,'datestr') & isfield(pdd,'timestr'))
   s{k}=['Date/time: ' pdd(1).datestr '/' pdd(1).timestr]; k=k+1;
end;
if isfield(pdd,'filename')
   s{k}=['File     : ' pdd(1).filename]; k=k+1;
end;   
if isfield(pdd,'fileid')
   s{k}=['File ID  : ' pdd(1).fileid]; k=k+1;
end;
if isfield(pdd,'comment')
   s{k}=['Comments : ' pdd(1).comment]; k=k+1;
end;
if isfield(pdd,'wavelength')
   s{k}=['Wavel.   : ' num2pstr(pdd(1).wavelength) 'm']; k=k+1;
end;
if isfield(pdd,'versionsensor')
   s{k}=['Sens. ver: ' num2pstr(pdd(1).versionsensor) ]; k=k+1;
end;


for j=1:n
  s{k}=['----------------------------------------']; k=k+1;
  s{k}=['Plane    : ' int2str(j)]; k=k+1;
  s{k}=['----------------------------------------']; k=k+1;
  
  if isfield(pdd,'z')  
     s{k}=['z        : ' num2str(pdd(j).z) ' [m]']; k=k+1;
  end;
  if isfield(pdd,'px')  
     s{k}=['pixels x : ' int2str(pdd(j).px)]; k=k+1;
  end;
  if isfield(pdd,'py')  
     s{k}=['pixels y : ' int2str(pdd(j).py)]; k=k+1;  
  end;
  if isfield(pdd,'focus')  
     s{k}=['focus    : ' int2str(pdd(j).focus)]; k=k+1;  
  end;
  if isfield(pdd,'rangexy')  
    s{k}=['xy range : ' num2str(pdd(j).rangexy) ' [m]']; k=k+1;
  end;
  if isfield(pdd,'windowsize')  
    s{k}=['win. size: ' num2str(pdd(j).windowsize) ]; k=k+1;
  end;
  if isfield(pdd,'posxy')  
     s{k}=['pos. xy  : ' num2str(pdd(j).posxy) ' [m]']; k=k+1;
  end;
  if isfield(pdd,'amp')    
     s{k}=['amplif.  : ' num2str(pdd(j).amp) ' [dB]']; k=k+1;  
  end  
  if isfield(pdd,'average')    
     s{k}=['averaging: ' int2str(pdd(j).average)]; k=k+1;    
  end;
  if isfield(pdd,'offset')      
     s{k}=['offset   : ' num2str(pdd(j).offset)]; k=k+1;      
  end;
  if isfield(pdd,'power')    
     s{k}=['power    : ' num2str(pdd(j).power)  ' [W]']; k=k+1;        
  end;
end; %for


if nargout<1 %If no output parameter is specified write to stdout
   for k=1:max(size(s))
       disp(s{k});
   end;
else
    cs=s;
end; %if



