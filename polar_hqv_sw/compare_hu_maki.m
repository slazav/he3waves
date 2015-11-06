function hqv()

  R=8;
  Lx=R+8;
  Ly=8;

  Nd=5; % points per xi

  BC=3; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0

  res = calc_text(R, Nd, Lx,Ly);
  [xxx,yyy]=meshgrid(res.xx, res.yy);
  A1 = hu_maki_func(xxx, yyy, res.R, 9.7, 0);
  A2 = hu_maki_func(xxx, yyy, res.R, 2.7, 1);

  t='Calculated texture vs Hu-Maki model'
  find_figure(t); clf; hold on; title(t);

  subplot(3,1,1); hold on; title('Calculated alpha_d');
    plot_angle(res.xx, res.yy, res.A);

  subplot(3,1,2); hold on; title('Hu-Maki F1');
    plot_angle(res.xx, res.yy, A1');

  subplot(3,1,3); hold on; title('Hu-Maki F2');
    plot_angle(res.xx, res.yy, A2');


end

function plot_angle(x,y,A)
  % plot data with a propper cut
  c=(length(y)+1)/2; %y central point
  A1 = A(:,c:end)';
  A2 = A(:,1:c)'; A2(find(A2<pi/2))=pi;
  surface(x, y(c:end), A1, 'edgecolor', 'none');
  surface(x, y(1:c), A2, 'edgecolor', 'none');
  %view(135,60);
  fixaxes();
  xlabel('x')
  ylabel('y')
  zlabel('alpha_d')
end
