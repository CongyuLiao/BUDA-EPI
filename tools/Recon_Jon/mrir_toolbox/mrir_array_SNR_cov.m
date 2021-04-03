function snr = mrir_array_SNR_cov(img_uncombined, Rn)
%MRIR_ARRAY_COMBINE_COV
%
% img_combine_cov = mrir_array_combine_cov(img_uncombined, Rn)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/26
% $Id: mrir_array_combine_cov.m,v 1.1 2008/01/27 00:27:29 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
 
  Ncol = mrir_ice_dimensions(img_uncombined, 'col');
  Nlin = mrir_ice_dimensions(img_uncombined, 'lin');
  Ncha = mrir_ice_dimensions(img_uncombined, 'cha');
  Npar = mrir_ice_dimensions(img_uncombined, 'par');
  
  % arbitrary constant (normalizes image to approx same level as RSS combo)
  alpha = sqrt(mean(diag(Rn)));
  
  dims = size(img_uncombined);
  dims(end+1:16) = 1;
  
  
  img_permute = permute(img_uncombined, [3, 1, 2, 4:ndims(img_uncombined)]);
  img_reshape = reshape(img_permute, Ncha, []);
  
  iRn = inv(Rn);
  
  snr_reshape = (img_reshape' * iRn) .* img_reshape.';

  snr = sqrt(abs( reshape(sum(snr_reshape, 2), [Ncol, Nlin, Ncha, dims(4:end)]) ));
  
  % recall that oversampling does not affect SNR map, but normalization
  % factor already applied to sensitivity map
  scale_factor = 1 / sqrt(2 * Ncol * Nlin * Npar);
  
  snr = scale_factor * snr;

  
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_cov.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
