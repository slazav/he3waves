function calc_ferm()

  % Calculation grid, square with Lx*Ly size, Nx*Ny points.

  Nx=85;
  Ny=85;
  Lx=400;
  Ly=400;

  kF = 0.18;   %
  m  = 1;   % particle mass
  xi = 150; % vortex size
  Delta=0.01;
%  Ef = kF/2/m;
%  xi=Ef/Delta;

  xx=linspace(-Lx/2, Lx/2, Nx);
  yy=linspace(-Ly/2, Ly/2, Ny);
  dx=Lx/(Nx-1);
  dy=Ly/(Ny-1);

  [xxx,yyy] = meshgrid(xx,yy);
  rrr = sqrt(xxx.^2 + yyy.^2);
  N=[Nx Ny];


  % differential operators:
  M = sparse(prod(N),prod(N));
  Dxx = M; Dyy = M; % d2/dx2, d2/dy2
  Dx  = M;  Dy = M; % d/dx, d/dy
  DD  = M; % diagonal part of the Delta operator

  % other space-dependent values
  D1 = Delta * (1-exp(-rrr/xi)) .* exp(1i*atan2(yyy,xxx));

  % fill matrices
  for ix=1:N(1)
    for iy=1:N(2)
      p0  = sub2ind(N, ix,iy);
      if ix>1    pxm = sub2ind(N, ix-1,iy); end
      if ix<Nx   pxp = sub2ind(N, ix+1,iy); end
      if iy>1    pym = sub2ind(N, ix,iy-1); end
      if iy<Ny   pyp = sub2ind(N, ix,iy+1); end

      % Dxx, Dyy (zero boundary cond)
      Dxx(p0,p0) = -2/dx^2;
      Dyy(p0,p0) = -2/dy^2;
      if (ix>1)  Dxx(p0,pxm) = 1/dx^2; end
      if (ix<Nx) Dxx(p0,pxp) = 1/dx^2; end
      if (iy>1)  Dyy(p0,pym) = 1/dy^2; end
      if (iy<Ny) Dyy(p0,pyp) = 1/dy^2; end

      % Dx,Dy (zero boundary cond)
      if (ix>1)  Dx(p0,pxm) = -1/(2*dx); end
      if (ix<Nx) Dx(p0,pxp) =  1/(2*dx); end
      if (iy>1)  Dy(p0,pym) = -1/(2*dy); end
      if (iy<Ny) Dy(p0,pyp) =  1/(2*dy); end

      % diagonal part of the Delta operator
      DD(p0,p0) = D1(iy,ix);
    end
  end

  % unit matrix
  I = diag(ones(Nx*Ny,1));

  % Kinetic energy
  E = ( -Dxx -Dyy - I*kF^2 )/(2*m);

  % Delta operator
  Do =  (-1i*Dx + Dy)/kF * DD;
  Dc =  (-1i*Dx - Dy)/kF * conj(DD);
  % 2x2 Hamiltonian:
  H = [ E Do ; Dc, -E];

  find_figure('a'); clf; hold on;
  %surface(xx,yy, abs(D1) )
  %spy(H)
  %return

  mgap=2.986e-4;
  [V,D] = eigs(H, 1, 0*mgap);
  
  nn=1
  U1 = reshape( V(1:Nx*Ny,nn), Nx,Ny);
  U2 = reshape( V(Nx*Ny+1:2*Nx*Ny,nn), Nx,Ny);
  D(nn,nn)

%  U1=U1./exp(1i*kF*rrr);
%  U2=U2./exp(1i*kF*rrr);

  m=max(max([abs(U1) abs(U2)]));
  clim=[-m m];
  h(1)=subplot(2,2,1);  hold on; title('real(U1)');
  surface(yy, xx, real(U1), 'EdgeColor','none')
  caxis(clim);
  h(2)=subplot(2,2,2);  hold on;
  surface(yy, xx, real(U2), 'EdgeColor','none')
  caxis(clim);
  h(3)=subplot(2,2,3);  hold on;
  surface(yy, xx, imag(U1), 'EdgeColor','none')
  caxis(clim);
  h(4)=subplot(2,2,4);  hold on;
  surface(yy, xx, imag(U2), 'EdgeColor','none')
  caxis(clim);

%  print /tmp/spwave2.png -noui -dpng

end