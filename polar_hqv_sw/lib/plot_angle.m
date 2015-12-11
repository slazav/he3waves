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
