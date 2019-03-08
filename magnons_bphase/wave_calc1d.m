% simplify wave_calc.m for 1D radial calculation

function [res sres] = wave_calc1d(cr,cz,bp,f0,fB,Imin,ctype, states, Nr)

  %  input:
  %    cr,cz,bp,f0,fB - perpendicular and parallel spin wave velocities, cm/s
  %    bp             - d(beta_N)/dr in the center
  %    f0, fB         - NMR frequncy, Leggett frequncym Hz
  %    Imin           - minimal coil current
  %    ctype:
  %       0 - analitycal formula
  %       1 - calculation
  %    states
  %    Nr,Nz - mesh size
  %
  %  output:
  %    res.rr  - 1d r-grid
  %    res.R   - cell radius, cm
  %    res.Nr  - Nr
  %    res.Urz - magnon potential, Hz above f0
  %
  %  output for each state:
  %    sres(ii).psi - normalized eigenfunction;
  %    sres(ii).df  - eigenfrequncy, Hz above f0
  %    sres(ii).Im  - integal Im = int(psi)^2/int(psi^2)
  %    sres(ii).Igr = sum(sum(2*pi*rrr.*Dr*u * dr*dz))^2 / Ie;
  %    sres(ii).Igz = sum(sum(2*pi*rrr.*Dz*u * dr*dz))^2 / Ie;


  R=0.2925; % cell radius, cm
  Mr = 1.6520; % Hmin coil paramters H = Mr(2z^2 - r^2)
  gyro = 20378;
  blp = bp * sqrt(5/2);
  e = (cz^2-cr^2)/(2*pi*f0) * blp;
  a = 2*(cz^2-cr^2)/cr^2 * blp;

  rr=linspace(0,R,Nr);
  dr=R/(Nr-1);

  res.rr  = rr;
  res.R = R;
  res.Nr = Nr;

  Ur2 = pi*fB^2/f0 * bp^2 - gyro*Mr*Imin; % quadratic terms in the energy
  res.Urz = Ur2*rr.^2;  % potential

  %%%%
  if ctype==0 % calculate analytically
    for nr=0:10

      wr=2*cr*sqrt(Ur2/(2*pi*f0));
      ar = cr/sqrt(pi*f0*wr);

      sres(nr+1).df  = wr*(2*nr+1)/(2*pi);
      sres(nr+1).psi = 1/(sqrt(pi)*ar) * exp(-rr.^2/(2*ar^2)) .* laguerrepoly(nr,(rr/ar).^2);

      if max(max(sres(nr+1).psi)) < -min(min(sres(nr+1).psi)); sres(nr+1).psi = -sres(nr+1).psi; end
      sres(nr+1).nr  = nr;
      sres(nr+1).nz  = 0;
      sres(nr+1).Im  = ar^2 * 8*pi^1.5;
      sres(nr+1).Igr  = 1/ar^2 * (2*nr+1);
    end
    return;
  end

  % we are building matrix for rr(1:Nr-1)
  % with Dirichlet boundary conditions at r(Nr)
  % and Neumann boundary condition at r(1).
  N=Nr-1;
  Drr = sparse(N,N);
  uDr = sparse(N,N); % 1/r d/dr
  Dr  = sparse(N,N); % d/dr
  U   = sparse(N,N);

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

  for ir=1:N % don't use Nr below!
      % Drr
      Drr(ir,ir) = -2/dr^2;
      if ir>1    Drr(ir,ir-1) = 1/dr^2; end
      if ir<N(1) Drr(ir,ir+1) = 1/dr^2; end
      if ir==1
        Drr(ir,ir+1) = 2/dr^2; % Neumann boundary condition
      end

      % uDr
      if ir>1
        uDr(ir,ir-1) = -1/(2*dr*rr(ir));
        if ir<N(1); uDr(ir,ir+1) =  1/(2*dr*rr(ir)); end
      else
        uDr(ir,ir) = -2/(dr^2);
        uDr(ir,ir+1) = 2/(dr^2);
      end

      % Dr
      if ir>1
        Dr(ir,ir-1) = -1/(2*dr);
        if ir<N(1); Dr(ir,ir+1) =  1/(2*dr); end
      end

      % U
      U(ir,ir)  = res.Urz(ir); % note +1 shift on the z axis
  end

  % non-modified equation
  A = - cr^2*(Drr + uDr)/(2*pi*f0) + U;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % solve eigenvalue problem
  [v,l] = eig(A);
  l=diag(l);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % classify states and calculate integrals
  for i=1:length(v(1,:)) % for every state
    % reshape result and add boundaries
    u=zeros(Nr,1);
    ur=zeros(Nr,1);
    u(1:Nr-1) = real(v(:,i));
    ur(1:Nr-1) = real(Dr*v(:,i));

    % find nz, nr - maybe not always accurate
    %number of zero crossings at r=0
    [~,im] = max(abs(u(1,:))); % max at r=0
    nr = length(find(u(1:end-1,im).*u(2:end,im)<0));

    % calculate integral for normalization
    Ie  = myint(rr, u.^2) * dr;

    % calculate other integrals
    sres(nr+1).Im  = (myint(rr, u)* dr)^2 / Ie;
    sres(nr+1).Igr = myint(rr, ur.^2)*dr / Ie;

    if max(max(u)) < -min(min(u)); u = -u; end
    sres(nr+1).df  = l(i)/(2*pi);
    sres(nr+1).psi = u/sqrt(Ie);
  end

end

%area integral
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
function res = myint(r, f)  % without dr factor
  dr=r(2)-r(1);
  fff = f(1:end-2)+f(3:end)-2*f(2:end-1);
  ff = f(1:end-2)-f(3:end);
  res = sum(r.*f')...
      - sum(r(2:end-1).* fff'/24)...
      + sum(dr*ff/12)...
      + f(1) * dr/8 ...
      + fff(1) * dr/128;
  res=res*2*pi;
end
