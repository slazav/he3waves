function hqv()

  % calculate lg in a box
  addpath lib

  find_figure('hqv wave 1');

  for R = [3.5 4.5 5.5] % soliton length
      for Nd=[4 6 8 10 12 14 16 18] % points per xi
        Lx=2*R;
        Ly=R;

        BC=0; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0
        res = calc_text(R, Nd, Lx,Ly);
        res=calc_wave(res,BC);

        i=1;
        %res.psi{i} = res.psi{i} .* exp(1i*res.A);
        A=sum(sum(res.psi{i}))^2 / sum(sum(res.psi{i}.^2));

        % print information
        ss=sprintf('%5.2f %4.1f %2d  %8.6f %8.6f %f\n', R, 0, Nd, ...
            real(res.en{1}), imag(res.en{1}), A);
        fo=fopen('lg_box.txt', 'a');
        fprintf(fo,ss);
        fclose(fo);
        fprintf(ss);

        % plot picture
        dx=(res.xx(end)-res.xx(1))*1.02;
        clf; hold on;
        title(ss);
        surface(res.xx,res.yy, real(res.psi{i}'), 'edgecolor', 'none');
        surface(res.xx+dx,res.yy, imag(res.psi{i}'), 'edgecolor', 'none');
        fixaxes();

        tt=sprintf('pics/pic_%05.2f_%02d.png', R, Nd);
        print('-dpng',tt);

      end
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
