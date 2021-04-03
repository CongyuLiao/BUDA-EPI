function snrmat = mrir_array_SNR_matrix(img, cov)
%MRIR_ARRAY_SNR_MATRIX
%
% snrmat = mrir_array_SNR_matrix(img, cov)

% the SNR matrix is used by mode mixing for SNR compression; the
% decomposition of the SNR matrix provides the transformation to the most
% SNR-efficient set of modes.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/09
% $Id: mrir_array_SNR_matrix.m,v 1.1 2008/11/09 22:37:05 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  W = mrir_array_whitening_operator(cov, 'svd');
  img_hat = mrir_array_transform(img, W);

  snrmat = sqrt(abs(mrir_array_stats_matrix(img_hat, 'cor')));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SNR_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
