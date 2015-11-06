function plot_lg()
  [R,m,N,lr,li] = textread('lg.txt',...
    '%f %f %f %f %f', 'commentstyle', 'shell');

  find_figure('lg'); clf; hold on;
  xlabel('grid density, pt/xi')
  ylabel('lg')

  RR=unique(R);
  mm=unique(m);

  R0=15

% for ir=3:length(RR)
%    R0=RR(ir);
    for im=1:length(mm)
      m0=mm(im);
      ii = find(R==R0 & m==m0);
      if length(ii)==0; continue; end
      x=N(ii);
      y=lr(ii)+1i*li(ii);
      [x,is]=sort(x); y=y(is);

      plot(x, real(y), 'r.-')
      plot(x, imag(y), 'b.-')
      text(x(end)+0.1, real(y(end)), sprintf('m=%.1f', m0));
      text(x(end)+0.1, imag(y(end)), sprintf('m=%.1f', m0));

      iif=find(x>9);
      if length(iif)>2
        ff1=fit(x(iif), real(y(iif)), 'a+b*exp(-x/c)', 'startpoint', [-1 1 10]);
        ff2=fit(x(iif), imag(y(iif)), 'a+b*exp(-x/c)', 'startpoint', [0 -1 10]);
        fprintf('%4.1f  %f %f\n', m0,ff1.a,ff2.a)
        xx=4:0.1:30;
        plot(xx, ff1(xx), 'k-')
        plot(xx, ff1.a, 'k--')
        plot(xx, ff2.a, 'k--')
      end
    end
%  end
end