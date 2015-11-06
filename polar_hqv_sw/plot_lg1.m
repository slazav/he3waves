function plot_lg()

  find_figure('lg1'); clf; hold on;

  %R/xi, real(lg), imag(lg)
  data1=[
    0.25 0 -1
    0.5  0 -1
    1 -0.399 -0.387
  1.8   -0.735 -0.212
    2 -0.763 -0.184
    3 -0.859 -0.104
    4 -0.902 -0.069
    5 -0.924 -0.049
    6 -0.939 -0.037
    7 -0.949 -0.029
    8 -0.955 -0.023
  ];
  data2=[
 2.0  -0.535 -0.090
 3.0  -0.589 -0.086
 3.5  -0.655 -0.088
 4.0  -0.741 -0.070
 4.5  -0.802 -0.065
 5.0  -0.854 -0.050
 5.5  -0.886 -0.047
 6.0  -0.913 -0.037
 7.0  -0.940 -0.029
 8.0  -0.953 -0.023
10.0  -0.966 -0.015
  ];


  xlabel('R/xi')
  ylabel('lg')

  x1=data1(:,1);
  y1=data1(:,2)+1i*data1(:,3);
  plot(x1, real(y1), 'r.-')
  plot(x1, imag(y1), 'b.-')

  x2=data2(:,1);
  y2=data2(:,2)+1i*data2(:,3);
  plot(x2, real(y2), 'r.-')
  plot(x2, imag(y2), 'b.-')

  xx=2:0.01:10;
  i1=find(x1>1.2);
  i2=find(x2>3);
  f1=fit(x1(i1), real(y1(i1)), 'a/x+b*exp(-x/c)', 'startpoint', [5 1 3]);
  f2=fit(x2(i2), real(y2(i2)), 'a/x+b*exp(-x/c)', 'startpoint', [5 1 3]);
  plot(xx,f1(xx), 'k-')
  plot(xx,f2(xx), 'k-')
  fprintf('%f + %f exp(-x/%f)\n', f1.a, f1.b, f1.c);
  fprintf('%f + %f exp(-x/%f)\n', f2.a, f2.b, f2.c);
%-0.939629 + 0.825742 * exp(-x/1.292818)
%-0.971340 + 3.151839 * exp(-x/1.525937)
%

end