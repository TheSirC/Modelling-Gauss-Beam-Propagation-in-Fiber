function pdduff=uffread(fullfilename)
%UFFREAD Read UFF file from disk
%   PDD=UFFREAD(FILENAME) reads the data from a file (D*.*) created by the 
%   "Laserscope" by PROMETEC GmbH representing the power density profile
%   (caustic) of a laser beam measured by the Laserscope UFF100 of 
%   PROMETEC GmbH (http://www.prometec.de) and returns a PDD stuct, with 
%   the following fields
%
%   name          (Prometec UFF100 Laserscope, string)
%   filename      (File ID, string)
%   comment       (Comments, string)        
%   px            (No of pixels along x-axis, integer)
%   py            (No of pixels along y-axis, integer)
%   focus         (Integer)
%   amp           (Amplification of the sensor signal, integer)
%   windowsize    (Window size, integer)
%   windowcenterx (x coordinate of window center, 48x1 double)
%   windowcentery (y coordinate of window center, 48x1 double)
%   versionsensor (Version sensor, integer)
%   versionuff    (Version UFF, integer)
%   power         (Laser power [W])
%   datestr       (Date, string)
%   timestr       (Time, string) 
%   x             (x coordinates, vector of px elements)
%   y             (y coordinates, vector of py elements)
%   ixy           (px*py matrix power density [W/m^2])
%
%   If FILENAME is ommitted, a dialog box for the user to select a filename
%   is displayed.
%  
%   EXAMPLE:
%   pdd=uffread('D1040180.005')
%   plotpdd(pdd)
%   disppdd(pdd)
%
%   See also mdfread

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

%Get filename and path if none were specified
if nargin~=1
       [filename,pathname]=uigetfile('D*.*','Select Prometec UFF file');
       fullfilename=[pathname filename];
    if fullfilename==0
       error('Laser Toolbox: UFFREAD: No file selected!','UFFREAD');
    end;
end;

if ~isempty(strfind(fullfilename,'.m'))
    warning(['Laser Toolbox: UFFREAD: MATLAB scripts are not valid ',...
        'PROMETEC UFF LASERSCOPE files. Press any key to continue.']);
    pause;
    
end;

[fid,message]=fopen(fullfilename,'r');

if fid<0
  error(['Laser Toolbox: uffread: File ' filename ' not found in ' pathname ]);
end;
disp(['Reading ' fullfilename ' ...']);


%Initialisation
pdd.name='Prometec UFF100 Laserscope';
pdd.filename=fullfilename;

fileinfo=dir(fullfilename);
pdd.datestr=datestr(fileinfo.datenum,1);
pdd.timestr=datestr(fileinfo.datenum,13);

gridbreite=80;
gridtiefe=40;

pdd.px=gridbreite+1;
pdd.py=gridtiefe+1;

%Reading data from fie
[prommat,count]=fread(fid,[gridbreite+1,gridtiefe+1],'uint8');
nullnum=fread(fid,48,'ubit1');
fokus=fread(fid,1,'uint8');
pdd.focus=fokus;
pdd.amp=fread(fid,1,'uint8'); %Verstaerkung
pdd.windowsize=fread(fid,1,'uint8'); %Fenstergroesse
pdd.windowcenterx=fread(fid,48,'ubit1'); %Fensterzentrum x
pdd.windowcentery=fread(fid,48,'ubit1'); %Fensterzentrum y
dummyexp=fread(fid,1,'uchar');
checksum=fread(fid,1,'long'); 
nochange=fread(fid,1,'uint8');
pdd.versionsensor=fread(fid,1,'int16');
pdd.versionuff=fread(fid,1,'int16');
fensterlaenge=fread(fid,1,'int16'); %Window length
fensterhoehe=fread(fid,1,'int16'); %Window heigth

if fokus==1
   gridsizeprom=fensterhoehe/1000;
else
   gridsizeprom=fensterhoehe;
end; %if fokus

pdd.power=fread(fid,1,'long'); %Power
pdd.comment=deblank(setstr(fread(fid,250,'uchar')')); %Comment

fclose(fid); %All data has been read

%Process some of the data 
deltabreite=gridsizeprom/(gridbreite);
deltatiefe=gridsizeprom/(gridtiefe);
breite=0:deltabreite:gridsizeprom;
tiefe=0:deltatiefe:gridsizeprom;
indexbreite=breite/deltabreite+1;
indextiefe=tiefe/deltatiefe+1;

pdd.x=(-gridbreite/2:gridbreite/2)*deltabreite*1e-3; 
pdd.y=(-gridtiefe/2:gridtiefe/2)*deltatiefe*1e-3; 

%Convert the 48 bits of nullnum (6 byte fl.point Real)
%Not yet used, by the way ...
enum=0;
for tel=1:8
	if nullnum(9-tel)==1
		enum=enum+2^(tel-1);
	end;
end;
s=nullnum(41);
fnum=0;
for bytetel=1:4
	for tel=1:8
		if nullnum(bytetel*8-tel+9)==1
			fnum=fnum+2^(tel-9+bytetel*8);
		end;
	end;
end;	
for tel=42:48
	if nullnum(90-tel)==1
		fnum=fnum+2^(tel-10);
	end;
end;
fnum=fnum/(2^39);
if enum==0
	nulldouble=0;
else
	nulldouble=(-1)^s*2^(enum-129)*(1+fnum);
end;	

%Calculate power density distribution in [W/m^2]
dx=deltatiefe*1E-3;
dy=deltabreite*1E-3;
m=min(min(prommat));
ixy=prommat-m;
pdd.ixy=(pdd.power*ixy/dx/dy/sum(sum(ixy)))';

if nargout==0
    disppdd(pdd);
else
    pdduff=pdd;
end;
