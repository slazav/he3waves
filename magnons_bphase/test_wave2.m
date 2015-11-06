function test_wave2()
  addpath /rota/programs/src/he3lib/matlab

  P=0;   % pressure
  ttc=0.1;    % temperature Tc
  cr = he3_cperp(ttc,P); % c_perp
  cz = he3_cpar(ttc,P);  % c_par
  fB = he3_nu_b(ttc,P);  % c_par
  bp = 2.5;
  f0 = 833000;
  Imin = 2;
  states = [00 02 04 10 20 30];
  plot_surf=0;

  Nr=30;
  Nz=30;

  [res1 sres1] = wave_calc(cr,cz,bp,f0,fB,Imin,0, states, Nr,Nz);
  [res2 sres2] = wave_calc(cr,cz,bp,f0,fB,Imin,1, states, Nr,Nz);
  [res3 sres3] = wave_calc(cr,cz,bp,f0,fB,Imin,2, states, Nr,Nz);

  figure; hold on;
  mm1=max(max(abs(sres1(1).psi)));
  mm2=max(max(abs(sres2(1).psi)));
  for i=1:6

    if plot_surf
      subplot(6,2,2*i-1); hold on;
      surface(res1.zz, res1.rr, sres1(i).psi);
      shading interp;
      xlim([-res1.Z res1.Z]/2);
      ylim([0 res1.R]);
      caxis([-mm1 mm1])

      subplot(6,2,2*i); hold on;
      surface(res2.zz, res2.rr, sres3(i).psi-sres2(i).psi);
      shading interp;
      xlim([-res2.Z res2.Z]/2);
      ylim([0 res2.R]);
      caxis([-mm2 mm2])
    else
      subplot(2,3,i); hold on;
      title(sprintf('%02d', states(i)));

      if length(sres1)>=i && prod(size(sres1(i).psi))
        plot(res1.zz, sres1(i).psi(1,:), 'k.-');
        plot(+res1.rr, sres1(i).psi(:,round(end/2)), 'k.-');
        plot(-res1.rr, sres1(i).psi(:,round(end/2)), 'k.-');
      end

      if length(sres2)>=i && prod(size(sres2(i).psi))
        plot(res2.zz, sres2(i).psi(1,:), 'r.--');
        plot(+res2.rr, sres2(i).psi(:,round(end/2)), 'r.--');
        plot(-res2.rr, sres2(i).psi(:,round(end/2)), 'r.--');
      end

      if length(sres3)>=i && prod(size(sres3(i).psi))
        plot(res3.zz, sres3(i).psi(1,:), 'g.--');
        plot(+res3.rr, sres3(i).psi(:,round(end/2)), 'g.--');
        plot(-res3.rr, sres3(i).psi(:,round(end/2)), 'g.--');
      end
      xlim([-res1.Z res1.Z]/2);
    end

    fprintf('%d%d -- %.2f %.2f %.2f  %.2f\n',...
      sres1(i).nr,sres1(i).nz, sres1(i).df, sres2(i).df, sres3(i).df, sres3(i).df/sres2(i).df);
  end

  kfr = (sres3(4).df - sres3(1).df)/(sres2(4).df - sres2(1).df);
  kfz = (sres3(2).df - sres3(1).df)/(sres2(2).df - sres2(1).df);
  fprintf('kfr= %.2f kfz=%.2f\n', kfr,kfz);

%  % plot result
%  figure;
%  surf(zz,rr,u);
end
