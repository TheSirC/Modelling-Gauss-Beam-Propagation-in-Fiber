%Laser Toolbox
%Version 0.1 beta (R2007b) September 20th, 2010
%
%Power density distribution.
%   gauss      - Gaussian power density profile.
%   templ      - Gauss-Laguerre mode.
%   temmn      - Gauss-Hermite mode.
%   tophat     - Top hat power density profile.
%   rectunif   - Rectangular uniform power density profile.
%   t2surfsrc  - Surface power density profile based on desired temperature profile.
%   mdfread    - Read MDF file from disk.
%   uffread    - Read UFF file from disk.
%   plotpdd    - Power density distribution plot.
%   disppdd    - Display power density distribution. 
%
%Laser beam.
%   iso11146   - Beam propagation ratios according to ISO 11146-1.
%   caustic    - Caustic plot.   
%   dispbeam   - Display beam characteristics.
%
%Temperature profile.
%   tpntsrc    - Temperature profile of point surface heat source.
%   tlinesrc   - Temperature profile of line heat source.
%   tsurfsrc   - Temperature profile of surface heat source.
%   plottemp   - Temperature plot.
%
%
%General.
%   overlap    - Overlap plot of pulses.
%   materials  - Sets and saves materials parameters.
%   momentxy   - First and second order moments.
%   num2pstr   - Convert numbers to a prefixed string.
%   hermite    - Hermite polynomial.
%   laguerre   - Laguerre polynomial.
%
%(c) 2008-2010 G.R.B.E. Römer, University of Twente.
%This Matlab Laser Toolbox is for private, non-commercial, single home 
%computer, or educational use only. The use of Matlab Laser Toolbox for 
%commercial purposes is strictly prohibited. Please read the detailed 
%license agreements in the "Matlab Laser Toolbox User Manual".
%
%
%