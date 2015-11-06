function hqv()

  R=50;
  Lx=70;
  Ly=70;

  Nd=1; % points per xi

  BC=3; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0
  M=1;

  res = calc_text(R, Nd, Lx,Ly, M);

  find_figure('hqv texture'); clf; hold on;

  subplot(3,1,1); hold on; title('Calculated alpha_d');
    plot_angle(res.xx,res.yy,res.A);

  subplot(3,1,2); hold on; title('(nabla alpha_d)^2');
    surface(res.xx,res.yy, res.xi*(res.Ax.^2+res.Ay.^2)', 'edgecolor', 'none');
    fixaxes();

  subplot(3,1,3); hold on; title('sin^2 alpha_d');
    surface(res.xx,res.yy, sin(res.A').^2, 'edgecolor', 'none');
    fixaxes();

  res=calc_wave(res,BC);

  find_figure('hqv wave'); clf; hold on;

  N=length(res.en);
  N=3;
  dx=(res.xx(end)-res.xx(1))*1.02;
  for i=1:N
    %res.psi{i} = res.psi{i} .* exp(1i*res.A);

    A=sum(sum(res.psi{i}))^2 / sum(sum(res.psi{i}.^2));
    ang=angle(res.psi{i}');

    if i==1;
      ii=find(ang<pi/2);
      ang(ii) = ang(ii)+2*pi;
      ang=ang-pi;
%    elseif i==2;
%      ii=find(ang(end/2:end,:)<0);
%      ang(ii) = ang(ii)-2*pi;
%    else
%      ii=find(ang<0.2*pi);
%      ang(ii) = ang(ii)+2*pi;
    end

    m=max(max(abs(res.psi{i})));
    subplot(N,1,i); hold on; title(sprintf('%.3f int:%f', res.en{i}, A));
      surface(res.xx,res.yy, abs(res.psi{i}')*pi/m, 'edgecolor', 'none');
      surface(res.xx+dx,res.yy, ang, 'edgecolor', 'none');
      fixaxes();
  end


end

function plot_angle(x,y,A)
  % plot data with a propper cut
  c=(length(y)+1)/2; %y central point
  A1 = A(:, c:end);
  A2 = A(:, 1:c);
  A2(:,c) = 2*pi-A1(:,1);
  surface(x,y(c:end), A1', 'edgecolor', 'none');
  surface(x,y(1:c), A2', 'edgecolor', 'none');

%  surface(x,y, A', 'edgecolor', 'none');
  %view(135,60);
  fixaxes();
  xlabel('x')
  ylabel('y')
  zlabel('alpha_d')
end
