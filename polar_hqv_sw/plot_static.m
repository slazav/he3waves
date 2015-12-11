function hqv()
  % plot various textures

  R=6;
  Lx=10;
  Ly=10;
  Nd=8; % points per xi
  BC=1; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0

  res = calc_text(R, Nd, Lx,Ly, 0.5);

  find_figure('hqv texture'); clf; hold on;


  subplot(3,2,1); hold on; title('alpha_d');
    plot_angle(res.xx,res.yy,res.A);
    caxis([0 2*pi]);
%    colorbar('location', 'northoutside');
  subplot(3,2,2); hold on; title('-sin^2 alpha_d');
    surface(res.xx,res.yy, -sin(res.A').^2, 'edgecolor', 'none');
    caxis([-1 0]);
%    colorbar('location', 'northoutside');
    fixaxes();

  res = calc_text(2*R, Nd, Lx,Ly, 0.5);

  subplot(3,2,3); hold on; title('alpha_d');
    plot_angle(res.xx,res.yy,res.A);
    caxis([0 2*pi]);
  subplot(3,2,4); hold on; title('-sin^2 alpha_d');
    surface(res.xx,res.yy, -sin(res.A').^2, 'edgecolor', 'none');
    caxis([-1 0]);
    fixaxes();

  res = calc_text(R, Nd, Lx,Ly, 1);

  subplot(3,2,5); hold on; title('alpha_d');
    plot_angle(res.xx,res.yy,res.A);
    caxis([0 2*pi]);
  subplot(3,2,6); hold on; title('-sin^2 alpha_d');
    surface(res.xx,res.yy, -sin(res.A').^2, 'edgecolor', 'none');
    caxis([-1 0]);
    fixaxes();


end
