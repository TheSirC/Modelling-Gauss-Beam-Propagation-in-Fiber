function [ dX,dY ] = calcul_ellipse( Px,Py,pt )
%calcul_ellipse Calculate the semi-minor axis and the semi-major axis ellipse at the given point
%   pt is the double at which the  computation takes place

Lx = find(Px(pt,:) >= 1/(exp(1)^2)); % Searching for the value corresponding to the beam diameter
Lx1 = Lx(1,1); % Lower bound for the beam diameter
Lx2 = Lx(1,end); % Upper bound for the beam diameter
dX = Lx2 - Lx1; % Beam diameter in pixel

Ly = find(Py(pt,:) >= 1/(exp(1)^2)); % Searching for the value corresponding to the beam diameter
Ly1 = Ly(1,1); % Lower bound for the beam diameter
Ly2 = Ly(1,end); % Upper bound for the beam diameter
dY = Ly2 - Ly1; % Beam diameter in pixel

end

