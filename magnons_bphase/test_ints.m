function test_ints()
  addpath /rota/programs/src/he3lib/matlab

  P=0;   % pressure
  ttc=0.1;    % temperature Tc
  cr = he3_cperp(ttc,P) % c_perp
  cz = he3_cpar(ttc,P)  % c_par
  fB = he3_nu_b(ttc,P)  % c_par
  bp = 3;
  f0 = 833000;
  Imin = 0.5;
  states = [00 02 04  10 12 14  20 22 24];
  plot_surf=0;

  figure; hold on;
  cols='krgbmc';
  leg={};

  h(1)=subplot(2,2,1); hold on; title('relative error in E');
  h(2)=subplot(2,2,2); hold on; title('relative error in Im');
  h(3)=subplot(2,2,3); hold on; title('relative error in Igr');
  h(4)=subplot(2,2,4); hold on; title('relative error in Igz');

  [res0 sres0] = wave_calc(cr,cz,bp,f0,fB,Imin,0, states,30,60);
  [res1 sres1] = wave_calc(cr,cz,bp,f0,fB,Imin,1, states,30,60);
  [res2 sres2] = wave_calc(cr,cz,bp,f0,fB,Imin,2, states,30,60);

  a0=[sres0.df]; a1=[sres1.df]; a2=[sres2.df];
  plot(h(1), 1:length(states), (a1-a0)./a0, 'ro-');
  plot(h(1), 1:length(states), (a2-a0)./a0, 'bo-');

  a0=[sres0.Im]; a1=[sres1.Im]; a2=[sres2.Im];
  plot(h(2), 1:length(states), (a1-a0)./a0, 'ro-');
  plot(h(2), 1:length(states), (a2-a0)./a0, 'bo-');

  a0=[sres0.Igr]; a1=[sres1.Igr]; a2=[sres2.Igr];
  plot(h(3), 1:length(states), (a1-a0)./a0, 'ro-');
  plot(h(3), 1:length(states), (a2-a0)./a0, 'bo-');

  a0=[sres0.Igz]; a1=[sres1.Igz]; a2=[sres2.Igz];
  plot(h(4), 1:length(states), (a1-a0)./a0, 'ro-');
  plot(h(4), 1:length(states), (a2-a0)./a0, 'bo-');

end
