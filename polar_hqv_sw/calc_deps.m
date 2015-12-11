function hqv()

  find_figure('hqv wave 1');
  addpath lib

  for R = [5] % soliton length
    for m = [-2]  % area size parameter
      for Nd=[12] % points per xi
        for t=[1] % solver type


          Lx=R*2+m;
          Ly=R*1+m;

          BC=1; % 0:periodic, 1:F=0, 2:F'=0, 3:x-per,y-F=0
          res = calc_text(R, Nd, Lx,Ly, 1);

          if t==1
            res = calc_wave(res,BC,0);
          elseif t==2
            res = calc_wave(res,BC,1);
          else
            res = calc_wave_q(res,BC);
          end
          %res.psi = res.psi .* exp(1i*res.A);
          A=abs(sum(sum(res.psi))*res.dx*res.dy)^2 /...
              sum(sum(abs(res.psi).^2)*res.dx*res.dy);

          % print information
          ss=sprintf('%5.2f %4.1f %2d %1d  %8.6f %8.6f %f\n', R, m, Nd, t, ...
                      abs(res.en), angle(res.en), A);
          fo = fopen('/rota/Analysis/temp/sla5/lg_amp.txt', 'a');
          fprintf(fo,ss);
          fclose(fo);
          fprintf(ss);

          % plot picture
          dx=(res.xx(end)-res.xx(1))*1.02;
          clf; hold on;
          title(ss);
          m1=max(max(abs(res.psi)));

%          surface(res.xx,res.yy, abs(res.psi')/m1, 'edgecolor', 'none');
%          surface(res.xx+dx,res.yy, real(res.psi')/m1, 'edgecolor', 'none');

          surface(res.xx,res.yy, abs(res.psi')/m1);
          surface(res.xx+dx,res.yy, real(res.psi')/m1);
          caxis([-1 1])

%        surface(res.xx,res.yy, imag(0.5i*sin(2*res.A').*res.psi'), 'edgecolor', 'none');
%        surface(res.xx+dx,res.yy, real(res.psi')/m1, 'edgecolor', 'none');
          fixaxes();

          tt=sprintf('/rota/Analysis/temp/sla5/pics/pic_%05.2f_%04.1f_%02d.png', R, m, Nd);
          print('-dpng',tt);
        end
      end
    end
  end


end

