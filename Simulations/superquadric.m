function retval=superquadric(epsilon,a)
  n=50;
  etamax=pi/2;
  etamin=-pi/2;
  wmax=pi;
  wmin=-pi;
  deta=(etamax-etamin)/n;
  dw=(wmax-wmin)/n;
  [i,j] = meshgrid(1:n+1,1:n+1)
  eta = etamin + (i-1) * deta;
  w   = wmin + (j-1) * dw;
  x = a(1) .* sign(cos(eta)) .* abs(cos(eta)).^epsilon(1) .* sign(cos(w)) .* abs(cos(w)).^epsilon(1);
  y = a(2) .* sign(cos(eta)) .* abs(cos(eta)).^epsilon(2) .* sign(sin(w)) .* abs(sin(w)).^epsilon(2);
  z = a(3) .* sign(sin(eta)) .* abs(sin(eta)).^epsilon(3);

  mesh(x,y,z);
  endfunction;