function T=num2pstr(x)

%NUM2PSTR Convert numbers to a prefixed string.
%    T = NUM2PSTR(x) converts the scalar x into a string representation T
%    with about 4 digits and a prefix of the metric system.  This is useful for
%    labeling plots with the TITLE, XLABEL, YLABEL, and TEXT commands.
%
%    See also num2str

%Author:	 Dr.ir. G.R.B.E. Römer, g.r.b.e.romer@utwente.nl
%Copyrights: All rigths reserved. G.R.B.E. Römer, University of Twente
%Date:		 10-june-2010


T=num2str(x);

%Prefix, see also http://en.wikipedia.org/wiki/SI_prefix
prefix=['y' 'z' 'a' 'f' 'p' 'n' 'u' 'm' ' ' 'k' 'M' 'G' 'T' 'P' 'E' 'Z' 'Y'];
i=[-24:3:-3 0 3:3:24]; %index to prefix
  for n=1:length(i)
      f=floor(abs(x)/10^i(n));
        if f<100
           T=[num2str(abs(x)/10^i(n)) prefix(n)];
           if x<0
               T=['-' T];
           end;
           return;
        end;
    end;
    

