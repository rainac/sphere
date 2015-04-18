function f = plotPattern3D()
  plot3([0], [0], [0]);

  set(gca, 'defaultlinelinewidth', 3);
  set(gca, 'xlim', [-0.1, 1.1]);
  set(gca, 'ylim', [-0.1, 1.1]);
  set(gca, 'zlim', [-0.1, 1.1]);
  
  set(gca, 'xtick', [0, 1]);
  set(gca, 'ytick', [0, 1]);
  set(gca, 'ztick', [0, 1]);

  set(gca, 'fontsize', 16);
  set(gca, 'defaulttextfontsize', 16);

  for i1=[0 1]
    for i2=[0 1]
      for j1=[0 1]
        for j2=[0 1]
          for k1=[0 1]
            for k2=[0 1]
              if (i1*i2 + j1*j2 + k1*k2) == 0 && ...
                    (i1 < i2  ...
                     || i1 == i1 && j1 < j2 ...
                     || i1 == i1 && j1 == j2 && k1 < k2)
                h = line([i1 i2], [j1 j2], [k1 k2]);
                if (i1 + i2 + j1 + j2 + k1 + k2) == 1
                  if i2 == 1
                    set(h, 'color', 'red');
                  elseif j2 == 1
                    set(h, 'color', 'green');
                  elseif k2 == 1
                    set(h, 'color', 'blue');
                  end
                elseif (i1 + i2 + j1 + j2 + k1 + k2) == 2
                  set(h, 'color', 'magenta');
                elseif (i1 + i2 + j1 + j2 + k1 + k2) == 3
                  set(h, 'color', [1 0.9 0.1]);
                else
                end
              elseif (i1 < i2  ...
                     || i1 == i1 && j1 < j2 ...
                     || i1 == i1 && j1 == j2 && k1 < k2)
                h = line([i1 i2], [j1 j2], [k1 k2]);
                set(h, 'color', 'black');
                set(h, 'linewidth', 1);
                set(h, 'linestyle', ':');
              end
            end
          end
        end
      end
    end
  end

  for i1=[0 1]
    for j1=[0 1]
      for k1=[0 1]
        if (k1 == 1)
          text(i1, j1, k1 + 0.05, sprintf('(%d %d %d)', i1, j1, ...
                                          k1));
        elseif (k1 == 0 && j1 == 0 && i1 == 0)
          text(i1, j1 + 0.1, k1 + 0.05, sprintf('(%d %d %d)', i1, j1, ...
                                          k1));
        else
          text(i1, j1, k1 - 0.05, sprintf('(%d %d %d)', i1, j1, ...
                                          k1));
        end
      end
    end
  end

  camorbit(150, -20);

  xlabel('x');
  hy = ylabel('y');
  cpos = get(hy, 'position');
  set(hy, 'position', cpos + [0, 0, 0.05]);
  zlabel('z');
  
  print('-depsc2', 'img/pattern-3d.eps');
  print('-r150', '-dpng', 'img/pattern-3d.png');
