function res=calc_wave(res,BC, tr)

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
  GA2 = sparse(prod(N),prod(N));
  SD  = sparse(prod(N),prod(N));
  SS  = sparse(prod(N),prod(N));
  JJ  = sparse(prod(N),prod(N));

  %vortex positions
  yv  = (res.Ny+1)/2;
  xv1 = (res.Nx+1)/2 - res.ixv;
  xv2 = (res.Nx+1)/2 + res.ixv;

  % "cut potential"
  J=zeros(size(res.A));
%  J(1:xv1,yv) = 1;
%  J(xv2:end,yv) = 1;

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


      % if f=a sqrt(x) + bx, f-1=0
      % f' = (sqrt(2) f1 - 3/2 f0)/(sqrt(2)-1)dx
      % f'' = - (sqrt(2) f1 - f0)/(sqrt(2)-1)4dx^2
      if iy==yv && ix==xv1+1
        DD(i0,i0) = 1/(4*(sqrt(2)-1)*dx^2) - 2/dy^2;
        DD(i0,ixp) = -sqrt(2)/(4*(sqrt(2)-1)*dx^2);
        DD(i0,ixm) = 0;
        Dx(i0,i0) = - 1.5/((sqrt(2)-1)*dx);
        Dx(i0,ixp) = sqrt(2)/((sqrt(2)-1)*dx);
        Dx(i0,ixm) = 0;
      end

      if iy==yv && ix==xv2-1
        DD(i0,i0) = 1/(4*(sqrt(2)-1)*dx^2) - 2/dy^2;
        DD(i0,ixm) = -sqrt(2)/(4*(sqrt(2)-1)*dx^2);
        DD(i0,ixp) = 0;
        Dx(i0,i0) = - 1.5/((sqrt(2)-1)*dx);
        Dx(i0,ixm) = sqrt(2)/((sqrt(2)-1)*dx);
        Dx(i0,ixp) = 0;
      end

      if iy==yv+1 && (ix==xv1 || ix==xv2)
        DD(i0,i0) = 1/(4*(sqrt(2)-1)*dx^2) - 2/dx^2;
        DD(i0,iyp) = -sqrt(2)/(4*(sqrt(2)-1)*dy^2);
        DD(i0,iym) = 0;
        Dx(i0,i0) = - 1.5/((sqrt(2)-1)*dy);
        Dx(i0,iyp) = sqrt(2)/((sqrt(2)-1)*dy);
        Dx(i0,iym) = 0;
      end

      if iy==yv-1 && (ix==xv1 || ix==xv2)
        DD(i0,i0) = 1/(4*(sqrt(2)-1)*dx^2) - 2/dx^2;
        DD(i0,iym) = -sqrt(2)/(4*(sqrt(2)-1)*dy^2);
        DD(i0,iyp) = 0;
        Dx(i0,i0) = - 1.5/((sqrt(2)-1)*dy);
        Dx(i0,iym) = sqrt(2)/((sqrt(2)-1)*dy);
        Dx(i0,iyp) = 0;
      end



      % F'=0 BC near the vortex
%      if iy==yv && ix==xv1
%        DD(i0,i0) = - 1/dx^2 - 2/dy^2;
%        DD(i0,ixp) = 1/dx^2;
%        DD(i0,ixm) = 0;
%        Dx(i0,i0) = - 1/2/dx;
%        Dx(i0,ixp) = 1/2/dx;
%        Dx(i0,ixm) = 0;
%      end
%
%      if iy==yv && ix==xv2
%        DD(i0,i0) = - 1/dx^2 - 2/dy^2;
%        DD(i0,ixm) = 1/dx^2;
%        DD(i0,ixp) = 0;
%        Dx(i0,i0) = 1/2/dx;
%        Dx(i0,ixm) = -1/2/dx;
%        Dx(i0,ixp) = 0;
%      end

      % GAx, GAy, SS, SS
      GAx(i0,i0) = res.Ax(ix,iy);
      GAy(i0,i0) = res.Ay(ix,iy);
      GA2(i0,i0) = res.Ax(ix,iy)^2 + res.Ay(ix,iy)^2;
      SD(i0,i0)  = sin(2*res.A(ix,iy));
      SS(i0,i0)  = sin(res.A(ix,iy))^2;
      JJ(i0,i0)  = 1e6*J(ix,iy);
    end
  end


  if tr
    A = JJ+ res.xi^2*(DD + (GAx*GAx + GAy*GAy)) + SS;
  else
    A = JJ+ res.xi^2*(DD + 2i*(GAx*Dx+GAy*Dy)) + SS + 1*0.5i*SD;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % solve eigenvalue problem 
  [v,l] = eigs(A,6,12);
  l=diag(l);

  [~,is] = sort(real(l));
  v=v(:,is);
  l=l(is);


  res.psi = reshape(v(:,end),Nx,Ny);
  res.en  = l(end);

  if tr; res.psi = res.psi .* exp(-1i*res.A); end
%  if ~tr; res.psi = res.psi .* exp(1i*res.A); end

%  res.psi=sin(res.A).^2;

end

