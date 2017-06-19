function plottemp(temp)
%PLOTTEMP Temperature plot.
%   PLOTTEMP(TEMP) plots a contour plot and a parametric mesh of the 3D
%   temperature profile at the minimum z-corrdinate specified in the struct
%   TEMP, as well as corresponding cross sections in the XZ plane and the 
%   YZ plane.
%
%   EXAMPLE
%   pdd=rectunif(35,[1e-3 2e-3],32)
%   plotpdd(pdd)
%   load C45
%   T=tsurfsrc(C45,pdd,50e-3,Inf,[0 0.01 0.02])
%   plottemp(T)
%
%   See also plotpdd, caustic

%   G.R.B.E. Römer 11-nov-2008
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

n=max(size(temp));
N=4; %Number of countours in contour plot
dig='%0.2E'; %format of data string (see SPRINTF for details)

for j=1:n %For each time instance a separate plot
    
    if max(max(temp(j).Txyz(:,:,1)))==0
        warning(['Laser Toolbox: PLOTTEMP: (rise of) temperature profile @ t=' num2str(temp(j).t) ' [s] is zero. Plot is skipped.']);
        continue;
    end;
        
        
    figure;
    if n>1
        set(gcf,'Name',[temp(j).name ', t=' num2str(temp(j).t) '[s]'  ]);
    else
       set(gcf,'Name',[temp(j).name ]);
    end;
    
    %Contour plot at z(1)
    subplot(221); 
    [cs,h]=contourf(temp(j).x,temp(j).y,temp(j).Txyz(:,:,1),N);
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    s='Contour';
    if isfield(temp(j),'z')
        s=[s ' @ z=' num2str(temp(j).z(1)) '[m]'];
    end;
    if n>1
        s=['t=' num2str(temp(j).t) '[s], '  s ];
    end;
    
    title(s);
    axis square;
    
    %y cross section
    subplot(222); 
    if (isfield(temp(j),'z') && (length(temp(j).z)>1) )
        sy=round(length(temp(j).y)/2);
        tT(:,:)=squeeze(temp(j).Txyz(:,sy,:));
        [cs,h]=contour(temp(j).y,-abs(temp(j).z),tT','*-.');
        clear tT;
        clabel(cs);
        ylabel('z [m]');
        xlabel('y [m]');  
        title(['Isotherms y-cross section @ x=' num2str(temp(j).y(sy)) 'm']);
    else
        subplot(222); 
        [M,i]=max(max(temp(j).Txyz(:,:)'));
        i=round(length(temp(j).x)/2); %
        plot(temp(j).Txyz(:,i),temp(j).y,'-k.');
        xlabel('y [m]');
        ylabel('T [K]');
        title(['Cross section @ x=' num2str(temp(j).x(i)) 'm.']);
        a=axis; a(2)=M*1.1; axis(a);
        axis square;
    end;
    
    %x cross section
    subplot(223); 
    if (isfield(temp(j),'z') && (length(temp(j).z)>1) )
        [M,i]=max(max(temp(j).Txyz(:,:,1)));
        sx=round(length(temp(j).x)/2);
        tT(:,:)=squeeze(temp(j).Txyz(sx,:,:));
        [cs,h]=contour(temp(j).x,-abs(temp(j).z),tT');
        clear tT;
        clabel(cs);
        ylabel('z [m]');
        xlabel('x [m]');
        title(['Isotherms x-cross section @ y=' num2str(temp(j).x(sx)) 'm']);
    else
        [M,i]=max(max(temp(j).Txyz(:,:)));
        i=round(length(temp(j).y)/2); %
        plot(temp(j).x,temp(j).Txyz(i,:),'-k.');
        xlabel('x [m]');
        ylabel('T [K]');
        title(['Cross section @ y=' num2str(temp(j).y(i)) 'm.']);
        a=axis; a(4)=M*1.1; axis(a);
        axis square;
    end;
    
    %isometric at z(1)
    subplot(224); 
    mesh(temp(j).x,temp(j).y,temp(j).Txyz(:,:,1));
    [m,i]=min(min(temp(j).Txyz(:,:,1)'));
    [M,i]=max(max(temp(j).Txyz(:,:,1)'));    
    a=axis; axis([a(1:4) m M]);
    axis vis3d;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('T [K]');
    if isfield(temp(j),'z')
       title(['Temperature profile @ z=' num2str(temp(j).z(1)) 'm']);
    else
       title(['Temperature profile']);
    end;

end; %for j