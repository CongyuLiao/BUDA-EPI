function i_repl = mrir_image_replicate(i_orig, dim, R)
%MRIR_IMAGE_REPLICATE  replicate image along one dimension, simulating ghosting or aliasing
%
% i_repl = mrir_image_replicate(i_orig, dim, R)

% intended for use with "mrir_array_GRAPPA_artifact_map.m"

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/01
% $Id: mrir_image_replicate.m,v 1.1 2008/05/01 04:55:59 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims = size(i_orig);
  dims(end+1:16) = 1;


  k_orig = mrir_fDFT_freqencode(mrir_fDFT_phasencode(i_orig));

  theta = zeros(dims(1), dims(2));

  for ind = 1:R,
    theta(ind:R:end,:) = pi/R*ind;
  end;

  k_repl = k_orig .* repmat(exp(i*theta.'), [1 1 dims(3:end)]);

  i_repl = mrir_iDFT_phasencode(mrir_iDFT_freqencode(k_repl));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_replicate.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: