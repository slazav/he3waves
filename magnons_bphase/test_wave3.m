function test_wave3()
  addpath /home/sla/he3lib/lib/matlab

  P=0;   % pressure
  ttc=0.1;    % temperature Tc
  cr = he3_text_cperp(ttc,P); % c_perp
  cz = he3_text_cpar(ttc,P);  % c_par
  fB = he3_nu_b(ttc,P);  % c_par
  bp = 3;
  f0 = 833000;
  states = [00 02 10];
  plot_surf=0;
  Imin=[200 500 1000 1500 2000 3000 4000]/1000;

  Nr=20;
  Nz=20;

  function [fr,fz] = get_f(cr,cz,bp,f0,fB,Imin,t, states, Nr,Nz)
    for i=1:length(Imin)
      [res sres] = wave_calc(cr,cz,bp,f0,fB,Imin(i),t, states, Nr,Nz);
      if length(sres)>2; fr(i) = sres(3).df-sres(1).df; else fr(i)=NaN; end
      if length(sres)>1; fz(i) = sres(2).df-sres(1).df; else fz(i)=NaN; end
    end
  end

%  [fr1,fz1] = get_f(cr,cz,bp,f0,fB,Imin,0, states, Nr,Nz);
  [fr2,fz2] = get_f(cr,cz,bp,f0,fB,Imin,1, states, Nr,Nz);
  [fr3,fz3] = get_f(cr,cz,bp,f0,fB,Imin,2, states, Nr,Nz);

  figure; hold on;
%  plot(Imin, fr1, 'ro-');
%  plot(Imin, fz1, 'ro-');
  plot(Imin, fr2, 'bo-');
  plot(Imin, fz2, 'bo-');
  plot(Imin, fr3, 'go-');
  plot(Imin, fz3, 'go-');

%  plot(Imin, fr1a, 'r--');
%  plot(Imin, fz1a, 'r--');
%  plot(Imin, fr2a, 'b*-');
%  plot(Imin, fz2a, 'b*-');


end
