function r = plotHistlist(fname)
  d=load(fname);

  d(:,2) = floor(d(:,2) ./ 0.1) * 0.1;

  figure
  plot(d(:,2), d(:,4), '.', 'markersize', 2);
  title('Particle height over time');
  xlabel('time (s)');
  ylabel('y-coord (m)');
  print('-deps2', 'img/hit-y-pos-bw.eps');
  print('-deps2c', 'img/hit-y-pos.eps');
  print('-dpng', '-r300', 'img/hit-y-pos.png');
  
  figure
  plot(d(:,2), d(:,5), '.', 'markersize', 2);
  title('Particle x-velocity over time');
  xlabel('time (s)');
  ylabel('x-vel (m/s)');
  print('-deps2', '-r300', 'img/hit-x-vel-bw.eps');
  print('-deps2c', '-r300', 'img/hit-x-vel.eps');
  print('-dpng', '-r300', 'img/hit-x-vel.png');

  figure
  plot(d(:,2), sqrt(d(:,5).^2 + d(:,6).^2), '.', 'markersize', 2);
  title('Particle absolute velocity over time');
  xlabel('time (s)');
  ylabel('velocity (m/s)');
  print('-dpng', '-r300', 'img/hit-abs-vel.png');

