function plot_lg()
  [R,m,N,lr,li, A] = textread('lg_amp.txt',...
    '%f %f %f %f %f %f', 'commentstyle', 'shell');

  find_figure('amp'); clf; hold on;
  xlabel('R/xi')
  ylabel('|Int(Psi)|^2/Int(|Psi|^2)')

%  [R,is]=sort(R); A=A(is);

  ii=find(R>7);
  p=polyfit(R(ii), A(ii),1);

  xx=[0 16];
  plot(R,A, 'r.-')
  plot(xx, polyval(p,xx), 'k-')
  plot(xx, 2*xx, 'k--')
  p(1)
-p(2)/p(1)
end