function img_full = mrir_array_SENSE(img_redu, sensitivity, covmtx, R)
%MRIR_ARRAY_SENSE
%
% img_full = mrir_array_SENSE_2d(img_redu, sensitivity, covmtx, R)
    
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/28
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtxinv = inv(covmtx);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  img_full = repmat(NaN, [Ncol, Nlin]);
  
  for jj = 1:Nlin/R,

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);
    
      for ii = 1:Ncol,
	% encoding matrix, E: Ncha x R
	E = permute(reshape( sensitivity(ii, aliased_pixels_ind, :), R, Ncha ), [2 1]);
	% noise-weighted pseudoinverse: Nalias x Nchan
	L = inv( (E' * covmtxinv * E) ) * E' * covmtxinv;
	
	% observed reduced-FOV image intensities: Nchan x 1
	s = squeeze(img_redu(ii, jj, :));
	
	% SENSE estimate: Nalias x 1
	m_hat = L * s;
	
	if ( any(~isfinite(m_hat)) ),
	  error('infinite');
	end;
	
	% stuff result into "img_full" matrix based on indices of aliased pixels
	img_full(ii, aliased_pixels_ind) = reshape(m_hat, [R, 1]);
      
      end;
      
  end;

  
  return;

  

  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
