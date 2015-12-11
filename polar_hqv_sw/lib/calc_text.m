function res = calc_text(R, Nd, Lx,Ly, M)
  % calculate static texture - pair of 2 vortixes with winding number M

  if nargin<5; M=0.5; end

  %  input:
  %    xi    -- xi_D parameter
  %    R     -- distance between HQV in xi
  %    Nd    -- number of points per xi
  %    Lx,Ly -- size of the calculation area (total) in xi

  %  calculation grid
  xi=1;
  R =R*xi/2;
  Lx=Lx*xi/2;
  Ly=Ly*xi/2;
  Nx=round(Nd*Lx/xi);
  Ny=round(Nd*Ly/xi);
  xx=linspace(0,Lx,Nx);
  yy=linspace(0,Ly,Ny);
  dx=Lx/(Nx-1);
  dy=Ly/(Ny-1);

  % vortex position
  ixv = find(xx>=R,1);
  if length(ixv)==0 ixv=length(xx)+1; end


  % differential operator:
  N=[Nx Ny];
  Dxx = sparse(prod(N),prod(N));
  Dyy = sparse(prod(N),prod(N));
  X  = zeros(prod(N),1);
  B  = zeros(prod(N),1);
  Bx = zeros(Nx,1);

  for ix=1:Nx
    for iy=1:Ny

      i0  = sub2ind(N, ix,iy);

      if iy==1 && ix<=ixv;
        B(i0) = M*pi;
        Bx(ix) = M*pi;
      end

      if ix>1  ixm = sub2ind(N, ix-1,iy); end
      if ix<Nx ixp = sub2ind(N, ix+1,iy); end
      if iy>1  iym = sub2ind(N, ix,iy-1); end
      if iy<Ny iyp = sub2ind(N, ix,iy+1); end

      Dxx(i0,i0) = -2/dx^2;
      Dyy(i0,i0) = -2/dy^2;
      if ix>1  Dxx(i0,ixm) = 1/dx^2; end
      if iy>1  Dyy(i0,iym) = 1/dy^2; end
      if ix<Nx Dxx(i0,ixp) = 1/dx^2; end
      if iy<Ny Dyy(i0,iyp) = 1/dy^2; end
      % boundary conditions
      if ix==1  Dxx(i0,ixp) = 2/dx^2; end % F'=0
      if ix==Nx Dxx(i0,ixm) = 1/dx^2; end % F=0
      if iy==Ny Dyy(i0,iym) = 1/dy^2; end % F=0
    end
  end


  F=@(x) (Dxx+Dyy)*x + B/dy^2 - sin(2*x)/2/xi^2;
  X = fsolve(F, X);
  X = reshape(X, N);

  if ixv>length(xx); res.R=R*2;
  else res.R   = xx(ixv)*2; end
  res.xi  = xi;
  res.dx  = dx;
  res.dy  = dy;
  res.ixv = ixv;

  % alpha_D angle
  res.Nx  = 2*Nx-1;
  res.Ny  = 2*Ny+1;
  res.xx  = [-xx(end:-1:2) xx];
  res.yy  = [-yy(end:-1:1)-dy 0 yy+dy];

  res.A  = horzcat(...
        vertcat( M*2*pi-X(end:-1:2,end:-1:1), M*2*pi-X(:,end:-1:1)),...
        vertcat( Bx(end:-1:2), Bx),...
        vertcat( X(end:-1:2,:), X));

%  res.A  = horzcat(...
%        vertcat( M*2*pi-X(end:-1:2,end:-1:1), M*2*pi-X(:,end:-1:1)),...
%        M*pi + zeros(res.Nx,1),...
%        vertcat( X(end:-1:2,:), X));

  % derivatives
  res.Ax = zeros(size(res.A));
  res.Ay = zeros(size(res.A));

  U=zeros(1,res.Ny);
  U(1:Ny)=2*pi*M;
  U(Ny+1)=pi*M;

  res.Ax(2:end-1,:) = res.A(3:end,:)-res.A(1:end-2,:);
  res.Ax(1,:) = res.A(2,:);
  res.Ax(end,:) = res.A(end-1,:);
  res.Ax=fix_der(M, res.Ax)/(2*dx);
  res.Ax(:,Ny+1)=0;

  res.Ay(:,2:end-1) = res.A(:,3:end)-res.A(:,1:end-2);
  res.Ay(:,1) = res.A(:,2);
  res.Ay(:,end) = res.A(:,end-1);
  res.Ay=fix_der(M, res.Ay)/(2*dy);



end

function G=fix_der(M,G)
  ii=find(G >  pi*M); G(ii) = G(ii)-2*pi*M;
  ii=find(G < -pi*M); G(ii) = G(ii)+2*pi*M;
end
