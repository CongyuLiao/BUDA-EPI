function [fig, ax] = mrir_array_GRAPPA_1d__plot_dat_acs(index_datlines, index_acslines, Ncol, Nlin)
%MRIR_ARRAY_GRAPPA_1D__PLOT_DAT_ACS
%
% [fig, ax] = mrir_array_GRAPPA_1d__plot_dat_acs(index_datlines, index_acslines, Ncol, Nlin)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/11
% $Id: mrir_array_GRAPPA_1d__plot_dat_acs.m,v 1.2 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  sample_vector = zeros(1,Nlin);
  sample_vector(index_datlines) = sample_vector(index_datlines) + 2;
  sample_vector(index_acslines) = sample_vector(index_acslines) + 1;
  sample_vector(Nlin/2+1)       = 4;


  fig = figure('name', mfilename);
  [ax, h1, h2] = plotyy(1,1,1,1, 'plot'); delete([h1, h2]);
  hold on;
  imagesc(sample_vector.');
  set(ax(1), 'XLim', [0.5, 1.5], 'XTick', linspace(0.5, 1.5, 9), 'XTickLabel', linspace(0, Ncol, 9));
  set(ax(2), 'XLim', [0.5, 1.5], 'XTick', []);

  set(ax(1), 'YLim', [1,Nlin], 'YTick', index_datlines(1:2:end), 'YTickLabel', index_datlines(1:2:end) - (Nlin/2) - 1);
  set(ax(2), 'YLim', [1,Nlin], 'YTick', index_datlines(1:2:end));

  ylabel(ax(1), 'ky');
  ylabel(ax(2), 'phase-encoded line number');


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_1d__plot_dat_acs.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
