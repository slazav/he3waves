function [u,v] = hu_maki_coord(x,y)
  % convert x,y to Hu-Maki coordinates u,v
  % x = cosh(u) cos(v)
  % y = sinh(u) sin(v)

  B = 1+x.^2+y.^2;
  D = B.^2 - 4*x.^2;

  v = sqrt(B/2 - sqrt(D)/2);
  u = abs(x./v);

  ii=find(v==0); % x=0 should be treated specially
  u(ii)=sqrt(1+y(ii).^2);

  ii=find(u<1); % rounding errors
  u(ii)=1;

  v=acos(v);
  u=acosh(u);
end
