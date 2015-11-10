function calc_ferm()

  % Volovik's equation for quasiclassical trajectory
  % in the orbital half-vortex in the polar phase.
  % According with the note of Nov 9, 2015.

  % trajectory parameters, b and theta
  b  = 0.3;
  th = pi/4;

  xi = 1; %  Delta/EF

  % one coordinate - distance along the trajectory
  n = 500;
  s = linspace(-40,40, n);
  ds=(s(end)-s(1))/(n-1);

  % rho and \bar\phi coordinates
  r   = sqrt(s.^2 + b.^2);
  cph = s./r;
  sph = b./r;

  % radial part of the gap
  D = 1-exp(-r/xi);

  % pauli matrices
  t1 = [0  1;  1  0];
  t2 = [0 -1;  1  0] * 1i;
  t3 = [1  0;  0 -1];

  % differential operator
  I = speye(n,n);
  E = sparse(2:n,1:n-1,1,n,n);
%  Ds = ( 2*sparse(2:n,1:n-1,1,n,n) ...
%       + 3*speye(n,n) ...
%       - 6*sparse(1:n-1,2:n,1,n,n) ...
%       + 1*sparse(1:n-2,3:n,1,n,n) )/(6*ds);

   Ds = (E-E')/(2*ds);
%  Ds = (E-I)/ds;
  Dss = (E+E'-2*I)/ds^2;

  U1 = D/2.*(cos(th)+cph);
  U2 = D/2.*(sin(th)+sph);

  % hamiltonian
  H = kron(t3, -1i*xi*Ds )...
    + kron(t1, diag( U1 ))...
    - kron(t2, diag( U2 ));

  opt.v0= [1-D 1-D]';
  nn=1;
  [V,E] = eigs(H,nn, -0.33 ,opt);
  E(nn,nn)
  f1=V(1:n,nn); f2=V(n+1:2*n,nn);

  find_figure('a'); clf; hold on;
  %spy(H)
  %return
  subplot(1,2,1); hold on; xlabel('s');
    plot(s, D, 'r');
    plot(s, U1, 'g');
    plot(s, U2, 'b');
    legend('Delta', 'U1', 'U2', 'location', 'north');
  subplot(1,2,2); hold on; xlabel('c');
    plot(s, abs(f1), 'r.-');
    plot(s, abs(f2), 'b.-');

%    plot(s, imag(f1), 'r.--');
%    plot(s, imag(f2), 'b.--');


%  print /tmp/swave2.png -noui -dpng

end