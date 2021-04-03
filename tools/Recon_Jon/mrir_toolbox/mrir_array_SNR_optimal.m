function snr = mrir_array_SNR_optimal(sensitivity, covmtx, varargin)
%MRIR_ARRAY_SNR_OPTIMAL
%
% snr = mrir_array_SNR_optimal(sensitivity, covmtx)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/08
% $Id: mrir_array_SNR_optimal.m,v 1.1 2007/06/01 21:48:34 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  FLAG__PREWHITEN = 0;

  
  %==--------------------------------------------------------------------==%
  
  if ( nargin > 3 ),
    intensities = varargin{1};
  else,
    intensities = sensitivity;
  end;
  
  
  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);
  Npar = mrir_ice_dimensions(sensitivity, 'par');

  if ( FLAG__PREWHITEN ),
    % whitening operator for pre-whitening data
    W = mrir_array_whitening_operator(covmtx, 'svd');

    sensitivity_whitened = reshape(sensitivity, [], size(W,1)) * W.';
    sensitivity = reshape(sensitivity_whitened, size(sensitivity));

    intensities_whitened = reshape(intensities, [], size(W,1)) * W.';
    intensities = reshape(intensities_whitened, size(intensities));
    
  else,
    covmtxinv = inv(covmtx);
  end;

  % preallocate
  snr = zeros(Ncol, Nlin);
  
  for ii = 1:Ncol,
    for jj = 1:Nlin,
      S = squeeze(sensitivity(ii, jj, :));
      I = squeeze(intensities(ii, jj, :));

      if ( FLAG__PREWHITEN ),
	noisepower_inv = S' * I;
      else,
	noisepower_inv = S' * covmtxinv * I;
      end;
      
      snr(ii, jj) = sqrt(abs( noisepower_inv ));
    end;
  end;

  % recall that oversampling does not affect SNR map, but normalization
  % factor already applied to sensitivity map
  scale_factor = 1 / sqrt(2 * Ncol * Nlin * Npar);
  
  snr = scale_factor * snr;

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SNR_optimal.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
