function calc_ferm()
  R=1;
  N=1000;
  m=1;

  rr=linspace(0, R, N);
  dr=R/(N-1);

  % differential operator: 
  N=[Nx Ny];
  Drr = sparse(2*N, 2*N);
  uDr = sparse(2*N, 2*N);
  A   = sparse(2*N, 2*N);

  X  = zeros(prod(N),1);
  B  = zeros(prod(N),1);
  Bx = zeros(Nx,1);

  for ir=1:N

    % fill two halfes of the matrix
    % note that all checks are done with ir, not i1
    for i1 = [ir N+ir]
      % Drr
      Drr(i1,i1) = -2/dr^2;
      if (ir>1 && ir<N)
        Drr(i1,i1-1) = 1/dx^2;
        Drr(i1,i1+1) = 1/dx^2;
      elseif (ir==1)
        Drr(i1,i1+1) = 2/dx^2;
      elseif (ir==N)
        Drr(i1,i1-1) = 1/dx^2;
      end

      % 1/r Dr
      if ir>1
        uDr(x0,xrm) = -1/(2*dr*rr(ir));
        if ir<N(1); uDr(x0,xrp) =  1/(2*dr*rr(ir)); end
      else
        uDr(x0,x0) = -2/(dr^2);
        uDr(x0,xrp) = 2/(dr^2);
      end



  end


end