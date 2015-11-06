function hqv()

  xi=1;
  R=5;
  Lx=R+8;
  Ly=8;

  Nd=10; % points per xi

  BC=3; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0

  res = calc_text(xi, R, Nd, Lx,Ly);

  find_figure('hqv texture'); clf; hold on;

  subplot(3,1,1); hold on; title('alpha_d');
    plot_angle(res);

  subplot(3,1,2); hold on; title('(nabla alpha_d)^2');
    surface(res.xx,res.yy, xi*(res.Ax.^2+res.Ay.^2)', 'edgecolor', 'none');
    fixaxes();

  subplot(3,1,3); hold on; title('sin^2 alpha_d');
    surface(res.xx,res.yy, sin(res.A').^2, 'edgecolor', 'none');
    fixaxes();

  res=calc_wave(res,BC);

  find_figure('hqv wave'); clf; hold on;

  N=length(res.en);
  dx=(res.xx(end)-res.xx(1))*1.02;
  for i=1:N
    %res.psi{i} = res.psi{i} .* exp(1i*res.A);

    A=sum(sum(res.psi{i}))^2 / sum(sum(res.psi{i}.^2));

    subplot(N/2,2,i); hold on; title(sprintf('%.3f int:%f', res.en{i}, A));
      surface(res.xx,res.yy, real(res.psi{i}'), 'edgecolor', 'none');
      surface(res.xx+dx,res.yy, imag(res.psi{i}'), 'edgecolor', 'none');
      fixaxes();
  end


end

function plot_angle(res)
  % plot data with a propper cuts
  c=(res.Ny+1)/2; %y central point
  A1 = res.A(:,c:end)';
  A2 = res.A(:,1:c)'; A2(find(A2<pi/2))=pi;
  surface(res.xx,res.yy(c:end), A1, 'edgecolor', 'none');
  surface(res.xx,res.yy(1:c), A2, 'edgecolor', 'none');
  %view(135,60);
  fixaxes();
  xlabel('x')
  ylabel('y')
  zlabel('alpha_d')
end
