function calc_ferm()

  % Volovik's equation for quasiclassical trajectory
  % in the orbital half-vortex in the polar phase.
  % According with the note of Nov 9, 2015.


  find_figure('a'); clf; hold on;
  %spy(H)
  %return
  res=calc_traj(0.5,0.8);
  res.x
  subplot(1,2,1); hold on; xlabel('s');
    plot(res.s, res.D, 'r');
    plot(res.s, res.U1, 'g');
    plot(res.s, res.U2, 'b');
    legend('Delta', 'U1', 'U2', 'location', 'north');
  subplot(1,2,2); hold on; xlabel('c');
    plot(s, abs(res.f1), 'r.-');
    plot(s, abs(res.f2), 'b.-');


  find_figure('c'); clf; hold on;

  if 1
    b=[-2 1 0.5 0 0.5 1 2];
    th=0:0.1:2*pi;
    for i=1:length(b)
      for j=1:length(th)
        res=calc_traj(b(i),th(j));
        e(j)=res.e;
        x(j)=res.x;
        fprintf('%f %f %f %i\n', b(i), th(j), e(j), res.x);
      end
      ii=find(x==0)
      plot(th, e, 'r');
      plot(th(ii), e(ii), 'b');
    end
  else
    b=-3:0.1:1.4;
    th=pi/2;
    for i=1:length(b)
      res=calc_traj(b(i),th);
      e(i)=abs(res.e);
      x(i)=res.x;
      fprintf('%f %f %f %i\n', b(i), th, e(i), res.x);
    end
    ii=find(x==0)
    plot(b, e, 'r.-');
    plot(b(ii), e(ii), 'b.-');
    plot(b, 0.8*(1-exp(-abs(b)/xi)), 'k-');
    plot(b, 0.1751./abs(b), 'k-');
  end
 



%      plot(s, imag(f1), 'r.--');
%      plot(s, imag(f2), 'b.--');


  function res=calc_traj(b,th)
    % calculate for trajectory parameters b and theta

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
%    Ds = ( 2*sparse(2:n,1:n-1,1,n,n) ...
%           + 3*speye(n,n) ...
%           - 6*sparse(1:n-1,2:n,1,n,n) ...
%           + 1*sparse(1:n-2,3:n,1,n,n) )/(6*ds);

      Ds = (E-E')/(2*ds);
%      Ds = (E-I)/ds;
      Dss = (E+E'-2*I)/ds^2;

      U1 = D/2.*(cos(th)+cph);
      U2 = D/2.*(sin(th)+sph);

      % hamiltonian
      H = kron(t3, -1i*xi*Ds )...
        + kron(t1, diag( U1 ))...
        - kron(t2, diag( U2 ));

      opt.v0= [1-D 1-D]';
      nn=1;
      [V,E] = eigs(H,nn, 0 ,opt);

      res.s  = s;
      res.D  = D;
      res.U1 = U1;
      res.U2 = U2;
      res.f1 = V(1:n,nn);
      res.f2 = V(n+1:2*n,nn);
      res.e  = E(nn,nn);
      if sum(abs([res.f1([1 2 end-1 end]) res.f2([1 2 end-1 end])])) > 0.005;
        res.x=1; else res.x=0;
      end
  end

%  print /tmp/swave2.png -noui -dpng

end