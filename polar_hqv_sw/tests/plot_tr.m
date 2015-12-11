function plot_tr()

  find_figure('tr'); clf; hold on;
  xlabel('grid density, pt/xi')
  ylabel('lg')

  h(1)=subplot(1,2,1); hold on;
  h(2)=subplot(1,2,2); hold on;

  plot_tr1(h, 'tr.txt', 'r.-');
  plot_tr1(h, 'tr2.txt', 'b.-');
end

function plot_tr1(h, fname,c)
  [R,m,N,t, lr,li,A] = textread(fname,...
    '%f %f %f %d %f %f %f', 'commentstyle', 'shell');


  i1=find(t==1);
  i2=find(t==2);

  plot(h(1), N(i1), lr(i1), c)
  plot(h(1), N(i2), lr(i2), c)

  plot(h(2), N(i1), A(i1), c)
  plot(h(2), N(i2), A(i2), c)

end
