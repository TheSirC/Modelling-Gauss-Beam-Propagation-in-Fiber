function pddmdf=mdfread(fullfilename)
%MDFREAD Read MDF file from disk
%   PDD=MDFREAD(FILENAME) reads the data from a file (extension *.mdf), 
%   created by the "LaserDiagnoseSoftware" (version 2.81) by PRIMES GmbH
%   representing the power density profile(s) (caustic) of a laser beam 
%   measured by the FocusMonitor of PRIMES GmbH (http://www.primes.de) and
%   returns an array of PDD stucts (one for each plane), with the following
%   fields
%
%   name     (FocusMonitor, string)
%   filename (Filename, string)
%   fileid   (File ID, string)
%   comment  (Comments, string)        
%   px       (No of pixels along x-axis, integer)
%   py       (No of pixels along y-axis, integer)
%   rangexy  (Row vector with measuring range in x [m] and y [m] direction)
%   z        (z-position [m] along the axis of propagation)
%   posxy    (Row vector with x [m] and y [m] position of measuring range)
%   amp      (Amplification of the sensor signal [dB], double)
%   average  (Number of averages, integer)
%   offset   (Offset (noise level) of data, integer)
%   lambda   (Wavelegth laser [m])
%   power    (Laser power [W])
%   foclen   (Focal length [m])
%   datestr  (Date, string)
%   timestr  (Time, string) 
%   dataxy   (Data, px*py matrix, integers)
%   x        (x coordinates, vector of px elements)
%   y        (y coordinates, vector of py elements)
%   dataxy   (Data, px*py matrix, integers)
%   ixy      (px*py matrix power density [W/m^2])
%
%   If FILENAME is ommitted, a dialog box for the user to select a filename
%   is displayed.
%  
%   EXAMPLE:
%   pdd=mdfread('example.mdf')
%   plotpdd(pdd)
%   disppdd(pdd)
%
%   See also uffread

%   G.R.B.E. Römer 26-nov-2009
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

Pdefault=1; %Default power P [W] is not specified in the file

%Get filename and path if none were specified
if nargin~=1
       [filename,pathname]=uigetfile('*.mdf','Select a Primes MDF file');
       fullfilename=[pathname filename];
    if fullfilename==0
       error('Laser Toolbox: MDFREAD: No file selected!','MDFREAD');
    end;
end;

[fid,message]=fopen(fullfilename,'r');

if fid<0
  error(['Laser Toolbox: uffread: File ' filename ' not found in ' pathname ]);
end;
disp(['Reading ' fullfilename ' ...']);

fileinfo=dir(fullfilename);

%Read the data of n planes (assuming n>=1)
n=1;
while ~feof(fid)
%Read the header of the plane
k=1; comments=1;
pdd(n).name='Primes FocusMonitor';
pdd(n).filename=fullfilename;
pdd(n).fileid =fgetl(fid);  %File identifier
pdd(n).comment=[];
while comments %Read lines with comments
    tline = fgetl(fid);
      if tline(1)==';'
        m(n).comment=[pdd(n).comment ' ' tline(2:length(tline))];
        k=k+1;
      else
        comments=0;  
      end;
end;

pix            = str2num(tline); %Number of pixels: in x-direction in y-direction
pdd(n).px        = pix(1); %Number of pixels: in x-direction in y-direction
pdd(n).py        = pix(2); %Number of pixels: in x-direction in y-direction
pdd(n).rangexy   = str2num(fgetl(fid))*1e-3; %Size of the measuring range: dimension in x (m) dimension in y (m)
pdd(n).z         = str2num(fgetl(fid))*1e-3; %POSZ Position along the axis of beam propagation: z-position (m)
pdd(n).posxy     = str2num(fgetl(fid))*1e-3; %Transverse position of the center of the measuring range: x-pos y-pos (m)
pdd(n).amp       = str2num(fgetl(fid)); %Amplification of the signal: Amplification (dB)
pdd(n).average   = str2num(fgetl(fid)); %Number of averages: number
pdd(n).offset    = str2num(fgetl(fid)); %Offset - value display from the device: Offset - value

if pdd(n).fileid=='MDF101'
    pdd(n).wavelength = str2num(fgetl(fid))*1e-3; %Wavelength [m]
    pdd(n).power     = str2num(fgetl(fid)); %Power [W]
    pdd(n).foclen    = str2num(fgetl(fid))*1e-3; %Focal length [m]
    pdd(n).datestr   = fgetl(fid); %Date (string)
    if isempty(pdd(n).datestr)
        pdd(n).datestr=datestr(fileinfo.datenum,1);
    end;
    pdd(n).timestr   = fgetl(fid); %time (string)
    if isempty(pdd(n).timestr)
        pdd(n).timestr=datestr(fileinfo.datenum,13);
    end;    
else
    pdd(n).power     = Pdefault; %No power defiend so default to 1
    warning(['Laser Toolbox: mdfread: no power defined in plane ' int2str(n) ', defaulting to ' num2str(Pdefault) ' W.']);
end;

%Read the intensity data of this plane
data=fscanf(fid,'%i',pdd(n).px*pdd(n).py);
k=1;
for i=1:pdd(n).px
    for j=1:pdd(n).py
       pdd(n).dataxy(i,j)=data(k);
       k=k+1;
    end;
end;

pdd(n).dataxy=rot90(pdd(n).dataxy);

tline=fgetl(fid); %Read newline character

%calculate intensity [W/m^2] and determine x and y axis of plane
deltax=pdd(n).rangexy(1)/pdd(n).px; % delta x
deltay=pdd(n).rangexy(2)/pdd(n).py; %delta y
dA=deltax*deltay;

pdd(n).x=linspace(-pdd(n).rangexy(1)/2,pdd(n).rangexy(1)/2,pdd(n).px); % x axis
%pdd(n).x=pdd(n).x+pdd(n).posxy(1);
pdd(n).y=linspace(-pdd(n).rangexy(2)/2,pdd(n).rangexy(2)/2,pdd(n).py); % y axis
%pdd(n).y=pdd(n).y+pdd(n).posxy(2);
data=pdd(n).dataxy-pdd(n).offset;
pdd(n).ixy=data/sum(sum(data))/dA*pdd(n).power; %Calculate intensity profile [W/m^2]
pdd(n).ixy=rot90(pdd(n).ixy,-1);
n=n+1; %let's go to the next plane

end; %while ~eof

if (n-1)==1
    disp(['1 plane was read.']);
else
   disp([int2str(n-1) ' planes were read.']);
end;
fclose(fid);

if nargout==0
    disppdd(pdd);
else
    pddmdf=pdd;
end;

