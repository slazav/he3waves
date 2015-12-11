function res=calc_wave(res,BC)

  dx=res.dx;
  dy=res.dy;
  Nx=(res.Nx+1)/2;
  Ny=(res.Ny+1)/2;

  N=[Nx Ny];
  DD  = sparse(prod(N),prod(N));
  Dx  = sparse(prod(N),prod(N));
  Dy  = sparse(prod(N),prod(N));
  GAx = sparse(prod(N),prod(N));
  GAy = sparse(prod(N),prod(N));
  SD  = sparse(prod(N),prod(N));
  SS  = sparse(prod(N),prod(N));
  GA2 = sparse(prod(N),prod(N));
  Bx = zeros(Nx,1);

  for ix=1:Nx
    for iy=1:Ny
      i0  = sub2ind(N, ix,iy);

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

      if ix==1;  DD(i0,ixp) = 2/dx^2; end
      if ix==1;  Dx(i0,ixp) = 0; end
      if iy==1 && ix<=res.ixv;  DD(i0,iyp) = 2/dy^2; end
      if iy==1 && ix<=res.ixv;  Dy(i0,iyp) = 0; end

      if iy==1 && ix>res.ixv;  DD(i0,iyp) = 0; end % symmetric BC
      if iy==1 && ix>res.ixv;  Dy(i0,iyp) = 1/dy; end

      % GAx, GAy, SS, SS
      xa = ix+Nx-1;
      ya = iy+Ny-1;
      GAx(i0,i0) = res.Ax(xa,ya);
      GAy(i0,i0) = res.Ay(xa,ya);
      GA2(i0,i0) = res.Ax(xa,ya)^2 + res.Ay(xa,ya)^2;
      SD(i0,i0)  = sin(2*res.A(xa,ya));
      SS(i0,i0)  = sin(res.A(xa,ya))^2;
    end
  end

  tr=1;
  if tr
    A = res.xi^2*(DD + GA2) + SS;
  else
    A = res.xi^2*(DD + 2i*(GAx*Dx+GAy*Dy)) + SS + 0.5i*SD;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % solve eigenvalue problem 
  [v,l] = eigs(A,6,1);
  l=diag(l);

  N=1;
  [~,is] = sort(real(l));
  v=v(:,is);
  l=l(is);

  X = reshape(v(:,end),Nx,Ny);
  res.en  = l(end);

  res.psi  = horzcat(...
        vertcat( X(end:-1:2,end:-1:1), X(:,end:-1:1)),...
        vertcat( X(end:-1:2,2:end), X(:,2:end)));

  if tr; res.psi = res.psi .* exp(1i*res.A); end

end

