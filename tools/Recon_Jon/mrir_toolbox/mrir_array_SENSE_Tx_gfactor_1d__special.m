function gfactor = mrir_array_SENSE_Tx_gfactor_1d__special(excitation, transmit, covmtx, R)
%MRIR_ARRAY_SENSE_TX_GFACTOR_1D__SPECIAL
%
% g = mrir_array_SENSE_Tx_gfactor_1d__special(excitation_profile, Btransmit, noisecovmtx, R)

% based on:
%
%  Zhu, Y., "RF Power Deposition and 'g-factor' in Parallel Transmit",
%  Proc. Intl. Soc. Mag. Reson. Med. 14 (2006)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/12
% $Id: mrir_array_SENSE_gfactor_1d__special.m,v 1.1 2007/05/05 03:35:28 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(transmit, 1); % frequency encoded
  Nlin = size(transmit, 2); % phase encoded
  Ncha = size(transmit, 3);

  % approximate noise covariance matrix on transmit side with that measured
  % on receive side, assuming transmit coils were used to measure noise data
  Phi_inv = inv(covmtx);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for jj = 1:Nlin/R,

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);

    for ii = 1:Ncol,
      % encoding matrix, C: Nreplicates x Ncha   (Nreplicates == Nalias == R)
      C = reshape( transmit(ii, aliased_pixels_ind, :), R, Ncha );

      errmtx_inv = C * Phi_inv * C';

%      if ( cond(errmtx_inv) > 1/sqrt(eps) ),
%        error('the error covariance matrix is singular!');
%      end;
%      errmtx = inv(errmtx_inv + 2*eye(size(errmtx_inv)));
      errmtx = inv(errmtx_inv);

      mu = diag(permute(excitation(ii, aliased_pixels_ind), [2, 1]));

      gfactor(ii, aliased_pixels_ind) = ...
	  sqrt(diag(abs( [mu' * errmtx * mu] ./ [mu' * inv(diag(diag(errmtx_inv))) * mu] )));

      % ("abs" not strictly needed, but roundoff produces tiny imaginary part to discard)
      
    end;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__special.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: