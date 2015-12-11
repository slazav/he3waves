function calc_ferm()

  Nx=200;
  Ny=200;
  Lx=1;
  Ly=1;

  lF = 0.012;  % 1/kF
  xi  = 0.14;   %
  xiv = 0.4;    % vortex size

  l = 0; % state quantum number

  stype = 1;  %superfluid type
    % 0 -- s-wave superconductor, 1Q vortex
    % 1 -- A-phase with l||z, 1Q vortex
    % 2 -- polar phase, orbital 1/2-quantum vortex

  % does not work!
  nonosc = 0;

  xx=linspace(-Lx, Lx, Nx);
  yy=linspace(-Ly, Ly, Ny);
  dx=2*Lx/(Nx-1);
  dy=2*Ly/(Ny-1);

  [xxx,yyy] = meshgrid(xx,yy);
  rrr = sqrt(xxx.^2 + yyy.^2);
  phi = atan2(yyy,xxx);
  N=Nx*Ny;

  % diagonal matrices and differential operators
  Ix = speye(Nx,Nx);
  Iy = speye(Ny,Ny);
  I = speye(N,N);
  Ex = sparse(2:Nx,1:Nx-1,1,Nx,Nx);
  Ey = sparse(2:Ny,1:Ny-1,1,Ny,Ny);
  L = kron((Ex+Ex'-2*Ix)/dx^2, Iy)...
    + kron(Ix, (Ey+Ey'-2*Iy)/dy^2);
  Dx = kron((Ex-Ex')/(2*dx), Iy);
  Dy = kron(Ix, (Ey-Ey')/(2*dy));

  function DU = spdiag(U)
    N=prod(size(U));
    DU=sparse(1:N,1:N, reshape(U,N,1),N,N);
  end

  % Kinetic energy
  if nonosc
    % e-eF = -1/2m [ ipF(f/r + r nabla f \hat r) + nabla^2 f ]
    DrD = spdiag(xxx./rrr)*Dx + spdiag(yyy./rrr)*Dy;
    rri = spdiag(1./rrr);
    E = -( L + 1i/lF*(rri + 2*DrD) ) * xi*lF/2;
  else
    E = -( L + I/lF^2 ) * xi*lF/2;
  end


  Drad = (1-exp(-rrr/xiv)); % radial part of the gap
  
  % 2x2 Hamiltonian:
  if stype==0
    % s-wave superconductor, 1Q vortex
    D0 = Drad .* exp(1i*phi);
    D  = spdiag(D0);
    H = [ E, D
             conj(D), -E];
  elseif stype==1
    % A-phase with l||z, 1Q vortex
    D0 = Drad .* exp(1i*phi);
    D  = spdiag(D0);
    H = [ E, (-1i*Dx + Dy)*lF*D
             (-1i*Dx - Dy)*lF*conj(D), -E];
  elseif stype==2
    % polar phase, orbital 1/2-quantum vortex
    D0 = Drad .* exp(1i*phi/2);
    D  = spdiag(D0);
    mx = spdiag(cos(phi/2));
    my = spdiag(sin(phi/2));
    u = -1i*(Dx*mx + Dy*my)*lF;
    H = [ E, u*D
             u*conj(D), -E];
  else
    fprintf('Wrong stype parameter\n');
    return;
  end

  mgap = 2*0.3636 *lF/xiv;
  nn=1;
  [V,E] = eigs(H, nn, l*mgap);
  U1 = reshape( V(1:Nx*Ny,nn), Nx,Ny);
  U2 = reshape( V(Nx*Ny+1:2*Nx*Ny,nn), Nx,Ny);
  E(nn,nn)

%  U1 = U1* exp(-1i*rrr/lF);
%  U2 = U2* exp(-1i*rrr/lF);

  find_figure('a'); clf; hold on;

  m=max(max([abs(U1) abs(U2)]));
  clim=[-m m];
  h(1)=subplot(2,2,1);  hold on; title('real(U1)');
    surface(yy, xx, real(U1), 'EdgeColor','none')
    caxis(clim);
  h(2)=subplot(2,2,2);  hold on; title('real(U2)');
    surface(yy, xx, real(U2), 'EdgeColor','none')
    caxis(clim);

  if (0)
    h(3)=subplot(2,2,3);  hold on; title('imag(U1)');
      surface(yy, xx, imag(U1), 'EdgeColor','none')
      caxis(clim);
    h(4)=subplot(2,2,4);  hold on; title('imag(U2)');
      surface(yy, xx, imag(U2), 'EdgeColor','none')
      caxis(clim);
  else
    h(3)=subplot(2,2,3);  hold on; title('imag(U1)');
      plot(xx, Drad(:,end/2)*m, 'k-');
      plot(xx, real(U1(:,end/2)), 'r-');
      plot(xx, imag(U1(:,end/2)), 'b-');
      ylim(clim)
    h(4)=subplot(2,2,4);  hold on; title('imag(U2)');
      plot(xx, Drad(:,end/2)*m, 'k-');
      plot(xx, real(U2(:,end/2)), 'r-');
      plot(xx, imag(U2(:,end/2)), 'b-');
      ylim(clim)
  end

%  print /tmp/swave2.png -noui -dpng

end