function res=calc_wave(res,BC)

  dx=res.dx;
  dy=res.dy;
  Nx=res.Nx;
  Ny=res.Ny;

  N=[Nx Ny];
  DD  = sparse(prod(N),prod(N));
  Dx  = sparse(prod(N),prod(N));
  Dy  = sparse(prod(N),prod(N));
  GAx = sparse(prod(N),prod(N));
  GAy = sparse(prod(N),prod(N));
  SD  = sparse(prod(N),prod(N));
  SS  = sparse(prod(N),prod(N));

  for ix=1:Nx
    for iy=1:Ny
      i0  = sub2ind(N, ix,iy);

      if BC==0 % periodic BC
        if ix>1  ixm = sub2ind(N, ix-1,iy); else ixm = sub2ind(N, Nx,iy); end
        if ix<Nx ixp = sub2ind(N, ix+1,iy); else ixp = sub2ind(N,  1,iy); end
        if iy>1  iym = sub2ind(N, ix,iy-1); else iym = sub2ind(N, ix,Ny); end
        if iy<Ny iyp = sub2ind(N, ix,iy+1); else iyp = sub2ind(N, ix, 1); end
        % DD
        DD(i0,i0) = - 2/dx^2 - 2/dy^2;
        DD(i0,ixm) = 1/dx^2;
        DD(i0,ixp) = 1/dx^2;
        DD(i0,iym) = 1/dy^2;
        DD(i0,iyp) = 1/dy^2;
        % Dx, Dy
        Dx(i0,ixm) = -1/2/dx;
        Dx(i0,ixp) =  1/2/dx;
        Dy(i0,iym) = -1/2/dy;
        Dy(i0,iyp) =  1/2/dy;

      elseif BC==3 % x:periodic BC y:F=0
        if ix>1  ixm = sub2ind(N, ix-1,iy); else ixm = sub2ind(N, Nx,iy); end
        if ix<Nx ixp = sub2ind(N, ix+1,iy); else ixp = sub2ind(N,  1,iy); end
        if iy>1  iym = sub2ind(N, ix,iy-1); end
        if iy<Ny iyp = sub2ind(N, ix,iy+1); end
        % DD
        DD(i0,i0) = - 2/dx^2 - 2/dy^2;
        DD(i0,ixm) = 1/dx^2;
        DD(i0,ixp) = 1/dx^2;
        if iy>1  DD(i0,iym) = 1/dy^2; end
        if iy<Ny DD(i0,iyp) = 1/dy^2; end
        % Dx, Dy
        Dx(i0,ixm) = -1/2/dx;
        Dx(i0,ixp) =  1/2/dx;
        if iy>1  Dy(i0,iym) = -1/2/dy; end
        if iy<Ny Dy(i0,iyp) =  1/2/dy; end

      elseif BC==1 % F=0 BC
        if ix>1  ixm = sub2ind(N, ix-1,iy); end
        if ix<Nx ixp = sub2ind(N, ix+1,iy); end
        if iy>1  iym = sub2ind(N, ix,iy-1); end
        if iy<Ny iyp = sub2ind(N, ix,iy+1); end

        % DD
        DD(i0,i0) = - 2/dx^2 - 2/dy^2;
        if ix>1  DD(i0,ixm) = 1/dx^2; end
        if ix<Nx DD(i0,ixp) = 1/dx^2; end
        if iy>1  DD(i0,iym) = 1/dy^2; end
        if iy<Ny DD(i0,iyp) = 1/dy^2; end

        % Dx, Dy
        if ix>1  Dx(i0,ixm) = -1/2/dx; end
        if ix<Nx Dx(i0,ixp) =  1/2/dx; end
        if iy>1  Dy(i0,iym) = -1/2/dy; end
        if iy<Ny Dy(i0,iyp) =  1/2/dy; end

      elseif BC==2 % F'=0
        if ix>1  ixm = sub2ind(N, ix-1,iy); end
        if ix<Nx ixp = sub2ind(N, ix+1,iy); end
        if iy>1  iym = sub2ind(N, ix,iy-1); end
        if iy<Ny iyp = sub2ind(N, ix,iy+1); end

        % DD
        DD(i0,i0) = - 2/dx^2 - 2/dy^2;
        if ix>1  DD(i0,ixm) = 1/dx^2; end
        if ix<Nx DD(i0,ixp) = 1/dx^2; end
        if iy>1  DD(i0,iym) = 1/dy^2; end
        if iy<Ny DD(i0,iyp) = 1/dy^2; end

        if ix==1;  DD(i0,ixp) = 2/dx^2; end
        if ix==Nx; DD(i0,ixm) = 2/dx^2; end
        if iy==1;  DD(i0,iyp) = 2/dy^2; end
        if iy==Ny; DD(i0,iym) = 2/dy^2; end

        % Dx, Dy
        if ix>1 && ix<Nx; Dx(i0,ixm) = -1/2/dx; Dx(i0,ixp) = 1/2/dx; end
        if iy>1 && iy<Ny; Dy(i0,iym) = -1/2/dy; Dy(i0,iyp) = 1/2/dy; end
      end

      % GAx, GAy, SS, SS
      GAx(i0,i0) = res.Ax(ix,iy);
      GAy(i0,i0) = res.Ay(ix,iy);
      SD(i0,i0)  = sin(2*res.A(ix,iy));
      SS(i0,i0)  = sin(res.A(ix,iy))^2;
    end
  end

  A = -res.xi^2*(DD + 2i*(GAx*Dx+GAy*Dy)) - SS - 0.5i*SD;

%  A = -res.xi^2*(DD + GAx*GAx+GAy*GAy) - SS;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % solve eigenvalue problem 
  [v,l] = eigs(A,6,-1);
  l=diag(l);

  l
  [~,is] = sort(real(l));
  v=v(:,is);
  l=l(is);

  for i=1:length(l)
    res.psi{i} = reshape(v(:,i),Nx,Ny);
    res.en{i}  = l(i);
  end

end

