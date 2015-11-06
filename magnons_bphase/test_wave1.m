function test_wave1()
  addpath /rota/programs/src/he3lib/matlab

  P=0;   % pressure
  ttc=0.1;    % temperature Tc
  cr = he3_cperp(ttc,P) % c_perp
  cz = he3_cpar(ttc,P)  % c_par
  fB = he3_nu_b(ttc,P)  % c_par
  bp = 3;
  f0 = 833000;
  Imin = 3;
  states = [00];
  plot_surf=0;

  figure; hold on;
  cols='krgbmc';
  leg={};

  h(1)=subplot(2,2,1); hold on;
  h(2)=subplot(2,2,2); hold on;
  h(3)=subplot(2,2,3); hold on;
  h(4)=subplot(2,2,4); hold on;

  for ctype=0:3
    [res sres] = wave_calc(cr,cz,bp,f0,fB,Imin,ctype, states,30,60);

    for i=1:length(states)

      fr=sres(i).psi(:,round(end/2));
      fz=sres(i).psi(1,:);
      if ctype==0; fr0=fr; fz0=fz; end

      plot(h(1), res.zz, fz, [cols(ctype+1) '.-']);
      xlabel('z')

      plot(h(2), +res.rr, fr, [cols(ctype+1) '.-']);
      plot(h(2), -res.rr, fr, [cols(ctype+1) '.-']);
      fprintf('%02d %1d -- %.2f\n', states(i), ctype, sres(i).df);
      xlabel('r')

      plot(h(3), res.zz,  fz-fz0, [cols(ctype+1) '.-']);
      plot(h(4), +res.rr, fr-fr0, [cols(ctype+1) '.-']);
      plot(h(4), -res.rr, fr-fr0, [cols(ctype+1) '.-']);
      leg{end+1} = sprintf('type: %d, state: %d', ctype, states(i));
    end
  end
  legend(h(1), leg);


%    fprintf('%02d -- %.2f %.2f %.2f\n', states(i), sres1(i).df, sres2(i).df, sres3(i).df);
%  end
%  xlim([-0.1 0.1]);
%  ylim([40 60]);

end
