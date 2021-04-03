function varargout = mrir_array_combine_opt(img_uncombined, sensitivity, covmtx)
%MRIR_ARRAY_COMBINE_OPT
%
% img_n = mrir_array_combine_optimalSNR(img_uncombined, sensitivity, covmtx)
%
%
% example:
%
%  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
%  sensitivity = mrir_sensitivity_map(double(raw), 100);
%  img_uncombined = mrir_conventional(double(raw));
%
%  img_n = mrir_array_combine_optimalSNR(img_uncombined, sensitivity, covmtx);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/26
% $Id: mrir_array_combine_optimalSNR.m,v 1.1 2007/01/28 03:15:11 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Psi_inv = inv(covmtx);

  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);

  % arbitrary constant (normalizes image to approx same level as RSS combo)
  alpha = sqrt(mean(diag(covmtx)));

  dims = size(img_uncombined);

  dims(3) = 1;

  % preallocate
  img = zeros(dims);
  img_norm_uniform_noise = zeros(dims);
  img_norm_uniform_signal = zeros(dims);

  for ii = 1:Ncol,
    for jj = 1:Nlin,
      for slc = 1:mrir_ice_dimensions(img_uncombined, 'slc'),
	
	for rep = 1:mrir_ice_dimensions(img_uncombined, 'rep'),

	  %                        1,  2, 3, 4, 5, 6, 7, 8, 9,   0
	b = squeeze(sensitivity(ii, jj, :, 1, 1, 1, rep, 1, 1, slc));
	
	lambda_n = alpha / sqrt( b' * Psi_inv * b );
	lambda_s = alpha /     ( b' * Psi_inv * b );
	
	w   =            (Psi_inv * b);
	w_n = lambda_n * (Psi_inv * b);
	w_s = lambda_s * (Psi_inv * b);
	
	  
	  %                           1,  2, 3, 4, 5, 6,   7, 8, 9,   0
	  s = squeeze(img_uncombined(ii, jj, :, 1, 1, 1, rep, 1, 1, slc));
	  
	  img(ii, jj, :, 1, 1, 1, rep, 1, 1, slc)                     = w'   * s;
	  
	  img_norm_uniform_noise(ii, jj, :, 1, 1, 1, rep, 1, 1, slc)  = w_n' * s;
	  img_norm_uniform_signal(ii, jj, :, 1, 1, 1, rep, 1, 1, slc) = w_s' * s;
	  
	end;
      end;
      
      
    end;
  end;
  
  if ( nargout > 0 ),

    % in general, the resulting combined image will be complex-valued, but
    % in the case where the sensitivity maps are the same as input images
    % the images will be effectively real-valued up to machine precision.
    varargout{1} = img_norm_uniform_noise;
    varargout{2} = img_norm_uniform_signal;
    varargout{3} = img;

  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/EPI_3D/MATLAB/mrir_array_combine_optimalSNR.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: