function test_wave1d()
  addpath ~/PROG/he3lib/octave

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

  Nr=200;

  [res1 sres1] = wave_calc1d(cr,cz,bp,f0,fB,Imin,0, states, Nr);
  [res2 sres2] = wave_calc1d(cr,cz,bp,f0,fB,Imin,1, states, Nr);

  figure; hold on;
size(res1.rr)
size(sres1(1).psi)
size(sres2(1).psi)

  for nr = 0:10

%    if length(sres1)>=nr && prod(size(sres1(nr).psi))
      plot(+res1.rr, -4*nr+sres1(nr+1).psi, 'k.-');
      plot(-res1.rr, -4*nr+sres1(nr+1).psi, 'k.-');
%    end
%    if length(sres2)>=nr && prod(size(sres2(nr).psi))
      plot(+res2.rr, -4*nr+sres2(nr+1).psi, 'r.--');
      plot(-res2.rr, -4*nr+sres2(nr+1).psi, 'r.--');
%    end
    fprintf('%d -- %.2f %.2f\n',...
      nr, sres1(nr+1).df, sres2(nr+1).df);

  end

end
