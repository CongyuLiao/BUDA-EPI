function [img, noisecov] = mrir_array_modemix_calibrate(raw, noise, phase)
%MRIR_ARRAY_MODEMIX_CALIBRATE  
%
% [img, noisecov] = mrir_array_modemix_calibrate(raw, noise, phase)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/oct/26
% $Id: mrir_array_modemix_calibrate.m,v 1.1 2008/11/06 05:45:51 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  % subtract off phase errors on each channel from image and noise data


  ndims = size(noise);
  phase_rep = reshape(phase, [1 1 ndims(3)]);

  ndims(3) = 1;
  phase_mat = exp(-i * repmat(deg2rad(phase_rep), ndims));

  noisecov = mrir_array_stats_matrix(noise .* phase_mat, 'cov', 1);


  %==---------------------------

  idims = size(raw);

  idims(3) = 1;
  phase_mat = exp(-i * repmat(deg2rad(phase_rep), idims));

  img = mrir_conventional(raw .* phase_mat);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_modemix_calibrate.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
