function snr = mrir_array_SNR_individual(sensitivity, covmtx)
%MRIR_ARRAY_SNR_INDIVIDUAL
%
% snr = mrir_array_SNR_individual(sensitivity, covmtx)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/08
% $Id: mrir_array_SNR_individual.m,v 1.1 2007/03/23 23:45:26 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);
  Npar = size(sensitivity, 9);

  zcov = diag(covmtx);

  % preallocate
  snr = zeros(Ncol, Nlin, Ncha);

  for ll = 1:Ncha,
    snr(:, :, ll) = abs( sensitivity(:,:,ll) / sqrt(zcov(ll)) );
  end;

  % recall that oversampling does not affect SNR map, but normalization
  % factor already applied to sensitivity map
  scale_factor = 1 / sqrt(2 * Ncol * Nlin * Npar);
  
  snr = scale_factor * snr;

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/COIL_ARRAYS/MATLAB/mrir_array_SNR_individual.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
