function [res sres] = wave_calc(cr,cz,bp,f0,fB,Imin,ctype, states, Nr, Nz)

  %  input:
  %    cr,cz,bp,f0,fB - perpendicular and parallel spin wave velocities, cm/s
  %    bp             - d(beta_N)/dr in the center
  %    f0, fB         - NMR frequncy, Leggett frequncym Hz
  %    Imin           - minimal coil current
  %    ctype:
  %       0 - analitycal formula without nonhomogeneous terms
  %       1 - calculation without non-homogeneous terms
  %       2 - calculation with 1st order non-homogeneous terms 
  %       3 - calculation with 2nd order non-homogeneous terms
  %    states         - states to calculate, (00,01,10,11 etc.)
  %    Nr,Nz - mesh size
  %
  %  output:
  %    res.zz  - 1d z-grid
  %    res.rr  - 1d r-grid
  %    res.R   - cell radius, cm
  %    res.Z   - cell, length, cm
  %    res.Nr  - Nr
  %    res.Nz  - Nz
  %    res.Urz - magnon potential, Hz above f0
  %
  %  output for each state:
  %    sres(ii).psi - normalized eigenfunction;
  %    sres(ii).df  - eigenfrequncy, Hz above f0
  %    sres(ii).Im  - integal Im = int(psi)^2/int(psi^2)

  %    sres(ii).Igr = sum(sum(2*pi*rrr.*Dr*u * dr*dz))^2 / Ie;
  %    sres(ii).Igz = sum(sum(2*pi*rrr.*Dz*u * dr*dz))^2 / Ie;
  %    sres(ii).nr  = nr;
  %    sres(ii).nz  = nz;


  R=0.2925; % cell radius, cm
  Z=0.6; % cell length
  Mr = 1.6520; % Hmin coil paramters H = Mr(2z^2 - r^2)
  gyro = 20378;
  blp = bp * sqrt(5/2);
  e = (cz^2-cr^2)/(2*pi*f0) * blp;
  a = 2*(cz^2-cr^2)/cr^2 * blp;

  rr=linspace(0,R,Nr);
  zz=linspace(-Z/2,Z/2,Nz);
  dr=R/(Nr-1);
  dz=Z/(Nz-1);

  res.zz  = zz;
  res.rr  = rr;
  res.R = R;
  res.Z = Z;
  res.Nr = Nr;
  res.Nz = Nz;

  [zzz,rrr] = meshgrid(zz,rr);
  Ur2 = pi*fB^2/f0 * bp^2 - gyro*Mr*Imin; % quadratic terms in the energy
  Uz2 = 2*gyro*Mr*Imin;
  res.Urz = Ur2*rrr.^2 + Uz2*zzz.^2;  % potential
  mm0 = [1 0 1/2 0 3/8 0 5/16 0 35/128 0 63/256 0];

  %%%%
  if ctype==0 % calculate analytically
    for i=1:length(states)
      nr=floor(states(i)/10);
      nz=states(i)-10*nr;

      wr=2*cr*sqrt(Ur2/(2*pi*f0));
      wz=2*cz*sqrt(Uz2/(2*pi*f0));

      ar = cr/sqrt(pi*f0*wr);
      az = cz/sqrt(pi*f0*wz);
      sres(i).df  = (wr*(2*nr+1) + wz*(nz+0.5))/(2*pi);

      sres(i).psi = (pi^1.5 * ar^2 * az * 2^nz * factorial(nz))^(-0.5) *...
                     exp(-rrr.^2/(2*ar^2) - zzz.^2/(2*az^2)) .*...
                     laguerrepoly(nr,(rrr/ar).^2) .*...
                     hermitepoly(nz,zzz/az);

      if max(max(sres(i).psi)) < -min(min(sres(i).psi)); sres(i).psi = -sres(i).psi; end

      sres(i).nr  = nr;
      sres(i).nz  = nz;
      sres(i).Im  = ar^2*az * 8*pi^1.5  * mm0(nz+1);
      sres(i).Igz  = 1/az^2 * (nz+1/2);
      sres(i).Igr  = 1/ar^2 * (2*nr+1);
    end
    return;
  end

  % we are building matrix for rr(1:Nr-1), zz(2:Nz-1)
  % with Dirichlet boundary conditions on z(1), z(Nz), r(Nr) 
  % and Neumann boundary condition at r(1).
  N=[Nr-1 Nz-2];
  Drr = sparse(prod(N),prod(N));
  Dzz = sparse(prod(N),prod(N));
  rDrz= sparse(prod(N),prod(N));
  r2Dzz= sparse(prod(N),prod(N));
  Dr  = sparse(prod(N),prod(N)); % d/dr
  uDr = sparse(prod(N),prod(N)); % 1/r d/dr
  Dz  = sparse(prod(N),prod(N));
  Vir = sparse(prod(N),prod(N)); % 1/r
  Vr  = sparse(prod(N),prod(N)); % r
  U   = sparse(prod(N),prod(N));
  U1   = sparse(prod(N),prod(N));
  U2   = sparse(prod(N),prod(N));

  % Differential operators on the grid.
  % r>0:
  %   f = ar^2 + br + c
  %
  %   f(r)    = f0  = ar^2 + b r + c
  %   f(r-dr) = fm  = ar^2 - 2 a r dr + a dr^2 + br - b dr + c
  %   f(r+dr) = fp  = ar^2 + 2 a r dr + a dr^2 + br + b dr + c
  %
  %   fp-fm = 2 (2ar + b) dr
  %   fm+fp-2f0 = 2a dr^2
  %
  %   uDr: [1/r df/dr] = 2a + b/r = (fp-fm)/(2 r dr)
  %   Drr: [d2f/dr2] = 2a = (fm+fp-2f0)/dr^2
  %
  % r=0:
  %   fp=fm
  %
  %   f = f0 + ar^2
  %   fp = f0 + a dr^2
  %   [1/r df/dr] = [2f/dr2] = 2a = 2(fp-f0)/dr^2

  % ir=1                   ir=N(2)
  % |-2 +2 ... +1 -2 +1 ... +1 -2|  - Drr *dr^2
  %
  % | 0  0 ... -1  0 +1 ... -1  0|  - Dr *2dr
  %
  % | 0  0 ... -1  0 +1 ... -1 0| - Drz * 4drdz
  % | 0  0 ...  0  0  0 ...  0 0|
  % | 0  0 ... +1  0 -1 ... +1 0|

  for ir=1:N(1) % don't use Nr, Nz below!
    for iz=1:N(2)
      x0  = sub2ind(N, ir,iz);

      if ir>1    xrm = sub2ind(N, ir-1,iz); end
      if ir<N(1) xrp = sub2ind(N, ir+1,iz); end
      if iz>1    xzm = sub2ind(N, ir,iz-1); end
      if iz<N(2) xzp = sub2ind(N, ir,iz+1); end

      % Drr,Dzz
      Dzz(x0,x0) = -2/dz^2;
      Drr(x0,x0) = -2/dr^2;
      if iz>1    Dzz(x0,xzm) = 1/dz^2; end
      if iz<N(2) Dzz(x0,xzp) = 1/dz^2; end
      if ir>1    Drr(x0,xrm) = 1/dr^2; end
      if ir<N(1) Drr(x0,xrp) = 1/dr^2; end
      if ir==1
        Drr(x0,xrp) = 2/dr^2; % Neumann boundary condition
      end

      if iz<N(2) && ir>1 && ir<N(1) rDrz(x0,sub2ind(N, ir+1,iz+1)) = +rr(ir+1)/(4*dr*dz); end
      if iz<N(2) && ir>1            rDrz(x0,sub2ind(N, ir-1,iz+1)) = -rr(ir-1)/(4*dr*dz); end
      if iz>1 && ir>1 && ir<N(1)    rDrz(x0,sub2ind(N, ir+1,iz-1)) = -rr(ir+1)/(4*dr*dz); end
      if iz>1 && ir>1               rDrz(x0,sub2ind(N, ir-1,iz-1)) = +rr(ir-1)/(4*dr*dz); end

      % r2Dzz
      r2Dzz(x0,x0) = -2*rr(ir)^2/dz^2;
      if iz>1    r2Dzz(x0,xzm) = rr(ir)^2/dz^2; end
      if iz<N(2) r2Dzz(x0,xzp) = rr(ir)^2/dz^2; end

      % Dz
      if iz>1;    Dz(x0,xzm) = -1/(2*dz); end
      if iz<N(2); Dz(x0,xzp) =  1/(2*dz); end

      % Dr
      if ir>1
        Dr(x0,xrm) = -1/(2*dr);
        if ir<N(1); Dr(x0,xrp) =  1/(2*dr); end
      end

      % uDr
      if ir>1
        uDr(x0,xrm) = -1/(2*dr*rr(ir));
        if ir<N(1); uDr(x0,xrp) =  1/(2*dr*rr(ir)); end
      else
        uDr(x0,x0) = -2/(dr^2);
        uDr(x0,xrp) = 2/(dr^2);
      end

%      Dr(x0,x0) = -1/dr;
%      if ir<N(1); Dr(x0,xrp) =  1/dr; end

      % U
      U(x0,x0)  = res.Urz(ir,iz+1); % note +1 shift on the z axis


      U1(x0,x0) = Ur2*rr(ir)^2 + Uz2*(zz(iz+1) + a*rr(ir)^2/4)^2;
      U2(x0,x0) = cz^2 - cr^2 * a^2*rr(ir)^2/4;

      Vr(x0,x0) = rr(ir);
      if ir>1; Vir(x0,x0) = 1./rr(ir); end
    end
  end

  % non-modified equation
  A = - (cr^2*(Drr + uDr) + cz^2*Dzz)/(2*pi*f0) + U;

  % additional first order terms
  if ctype==2
    A = A - 2*e*(Dz + rDrz);
%    A = - (cr^2*(Drr + uDr + a*rDrz + a*Dz) + cz^2*Dzz)/(2*pi*f0) + U;
  end

  % additional first+second order terms
  if ctype==3
    A = A - e * (...
         2*(Dz + rDrz) ...
         + blp*(3*Vr*Dr + Vr*Vr*Drr) - blp * Vr*Vr*Dzz);
  end

  % equation after variable change
  if ctype==4
     A = - (cr^2*(Drr + uDr) + U2*Dzz)/(2*pi*f0) + U1;
  end

%Drr + 1/r Dr  ->  Drr + 1/r Dr - a rDrz - a Dz + (ar/2)^2 Dzz


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % solve eigenvalue problem
  [v,l] = eig(A);
  l=diag(l);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % classify states and calculate integrals
  for i=1:length(v(1,:)) % for every state
    % reshape result and add boundaries
    u=zeros(Nr,Nz);
    uz=zeros(Nr,Nz);
    ur=zeros(Nr,Nz);
    u(1:Nr-1,2:Nz-1)=reshape(real(v(:,i)),Nr-1,Nz-2);
    uz(1:Nr-1,2:Nz-1)=reshape(real(Dz*v(:,i)),Nr-1,Nz-2);
    ur(1:Nr-1,2:Nz-1)=reshape(real(Dr*v(:,i)),Nr-1,Nz-2);

    % find nz, nr - maybe not always accurate
    %number of zero crossings at r=0
    nz = length(find(u(1,1:end-1).*u(1,2:end)<0));
    [~,im] = max(abs(u(1,:))); % max at r=0
    nr = length(find(u(1:end-1,im).*u(2:end,im)<0));

    % we need only some of nr,nz
    if nr>=10 || nz>=10; continue; end
    ii=find(states == nr*10+nz);
    if length(ii)==0; continue; end

    % calculate integral for normalization
    Ie  = myint(rrr, u.^2) * dr*dz;

    % calculate other integrals
    sres(ii).Im  = (myint(rrr, u)* dr*dz)^2 / Ie;
    sres(ii).Igr = myint(rrr, ur.^2)* dr*dz / Ie;
    sres(ii).Igz = myint(rrr, uz.^2)* dr*dz / Ie;

    if max(max(u)) < -min(min(u)); u = -u; end
    sres(ii).df  = l(i)/(2*pi);
    sres(ii).psi = u/sqrt(Ie);
    sres(ii).nr  = nr;
    sres(ii).nz  = nz;

  end

end

%volume integral
  % r!=0:
  %   f = ar^2 + br + c
  %
  %   Vi/2pi = int(r f, r-dr/2,r+dr/2) = 
  %      [a(r+dr/2)^4/4 + b(r+dr/2)^3/3 + c(r+dr/2)^2/2] -
  %      [a(r-dr/2)^4/4 + b(r-dr/2)^3/3 + c(r-dr/2)^2/2] =
  %     f0 r dr - (fm+fp-2f0)/24 r dr  + (fp-fm)/12/dr^2
  %
  % r=0:
  %   f = f0 + ar^2
  %
  %   Vi/2pi = 2 int(r f, 0, dr/2) =
  %         [f0 dr^2/8 + (fm+fp-2f0) dr^2/128];
function res = myint(r, f)  % without dr*dz factor
  dr=r(2,1)-r(1,1);
  fff = f(1:end-2,:)+f(3:end,:)-2*f(2:end-1,:);
  ff = f(1:end-2,:)-f(3:end,:);
  res = sum(sum(r.*f))...
      - sum(sum(r(2:end-1,:).* fff/24))...
      + sum(sum(dr*ff/12))...
      + sum(f(1,:) * dr/8)...
      + sum(fff(1,:) * dr/128);
  res=res*2*pi;
end
