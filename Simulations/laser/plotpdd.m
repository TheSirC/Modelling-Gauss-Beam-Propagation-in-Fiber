function plotpdd(pdd)
%PLOTPDD Power density distribution plot.
%   PLOTPDD(PDD) plots a contour plot, two cross section as well as a 
%   parametric mesh of the power density profile(s) in the struct PDD; one
%   figure for each plane.
%
%   EXAMPLE
%   pdd=mdfread('example.mdf')
%   plotpdd(pdd)
%
%   See also plottemp, caustic

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente


if nargin~=1
    error('Laser Toolbox: PLOTPDD: too many or no input arguments.')
end;

%Initialisation
n=max(size(pdd));
N=4; %Number of countours in contour plot
dig='%0.2E'; %format of data string (see SPRINTF for details)
cmap='default'; %default colormap (see GRAPH3D for more maps)

if n>25
  ButtonName=questdlg(['This will create ' int2str(n) ' power density figures. Continue?'], ...
                      'Laser Toolbox: plotpdd', ...
                      'Yes','No','No');
   if strcmp(ButtonName,'No')
       return;
   end;
end; %if

for j=1:n
    Io=max(max(pdd(j).ixy)); %Peak intensity
    colormap(cmap);
%    pdd(j).ixy=rot90(pdd(j).ixy,-1);

   %Contour plot
    subplot(221); 
    
    [cs,h]=contourf(pdd(j).x,pdd(j).y,pdd(j).ixy,N);
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    title(['I_{peak}=' num2pstr(Io) 'W/m^2, P=' num2pstr(pdd(j).power) 'W.']);
    axis equal;

    %y cross section
    subplot(222); 
    [M,i]=max(max(pdd(j).ixy'));
    i=round(length(pdd(j).x)/2); %
    plot(pdd(j).ixy(:,i),pdd(j).y,'-k.');
    ylabel('y [m]');
    xlabel('I [W/m^2]');
    title(['Cross section @ x=' num2pstr(pdd(j).x(i)) 'm.']);
    a=axis; a(2)=Io*1.1; axis(a);
    axis square;
    
     %x cross section
    subplot(223); 
    [M,i]=max(max(pdd(j).ixy));
    i=round(length(pdd(j).y)/2); %
    plot(pdd(j).x,pdd(j).ixy(i,:),'-k.');
    xlabel('x [m]');
    ylabel('I [W/m^2]');
    title(['Cross section @ y=' num2pstr(pdd(j).y(i)) 'm.']);
    a=axis; a(4)=Io*1.1; axis(a);
    axis square;

    %isometric
    subplot(224); 
    mesh(pdd(j).x,pdd(j).y,pdd(j).ixy);
    [m,i]=min(min(pdd(j).ixy'));
    a=axis; axis([a(1:4) m M]);
    axis vis3d;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('I [W/m^2]');
    
    colorbar;
    s='';
    if isfield(pdd(j),'z') 
        s=['@ z=' num2str(pdd(j).z) 'm.'];
        title(s);
        set(gcf,'name',s);
    end;
    
    if isfield(pdd(j),'name')
       set(gcf,'Name',[pdd(j).name s]);
    end;

       
    if ((n>1) && (j<n))
        figure; %Create new figure if there 2 planes or more
    end;
    
end;


