function alpha = hu_maki_func(x,y, R, xi, N)
  % model texture used in Hu Maki paper
  % N = 0 or 1

  [u,v] = hu_maki_coord(2*x/R,2*y/R);
  if N==0 F=exp(-xi*(cosh(u)-1));
  else    F=exp(-xi*sinh(u)); end

  alpha = atan(F.*sin(v)./sinh(u));
  ii=find(y==0 & abs(x)<R/2);
  alpha(ii)=pi/2;

  ii=find(y==0 & abs(x)>=R/2);
  alpha(ii)=0;

  ii=find(y<0);
  alpha(ii)=pi-alpha(ii);
end
