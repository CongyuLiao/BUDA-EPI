function fig = mrir_array_GRAPPA_1d__plot_kernel(X1, Y1, X2, Y2, Ncol, Nlin, Nsrc_FE, Nsrc_PE, ii, jj)

  
  fig = 999;
  figure(fig);

  if ( ~any(ismember(allchild(0), fig)) ),
    axis; hold on; box on;
  end;
  
  
  cla;
  plot(X1(:), index_datlines(Y1(:)), 'ko', 'MarkerFaceColor', 'k');
  plot(X2(:), index_acslines(Y2(:)), 'wo', 'MarkerFaceColor', 'w');
  set(gca, 'XLim', [-Nsrc_FE,Ncol+Nsrc_FE], 'YLim', [-1+Nacs, Nlin-Nacs], 'Color', [1,1,1]*0.85);
  set(gca, 'YDir', 'normal');
  %	set(gca, 'YTick', linspace(0, Nlin, 33)+1);
  %	set(gca, 'YTickLabel', linspace(-Nlin/2, +Nlin/2, 33))
  xlabel('frequency-encoded');
  ylabel('phase-encoded');
  title(sprintf('ii = %d, jj = %d', ii, jj)); axis equal;
