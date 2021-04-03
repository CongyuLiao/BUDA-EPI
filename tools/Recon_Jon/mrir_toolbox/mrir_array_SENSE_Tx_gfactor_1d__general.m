function gfactor = mrir_array_SENSE_Tx_gfactor_1d__general(excitation, transmit, covmtx, R)
%MRIR_ARRAY_GFACTOR_1D__GENERAL
%
% g = mrir_array_SENSE_Tx_gfactor_1d__general(excitation_profile, Btransmit, noisecovmtx, R)
  
% based on:
%
%  Zhu, Y., "RF Power Deposition and 'g-factor' in Parallel Transmit",
%  Proc. Intl. Soc. Mag. Reson. Med. 14 (2006)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/12
% $Id: mrir_array_SENSE_gfactor_1d__general.m,v 1.2 2007/05/14 02:36:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(transmit, 1); % frequency encoded
  Nlin = size(transmit, 2); % phase encoded
  Ncha = size(transmit, 3);

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(covmtx, 'svd');

  transmit_whitened = reshape(transmit, [], Ncha) * W.';
  transmit = reshape(transmit_whitened, size(transmit));

  % image folding operator, F: ceil(Nlin/R) x Nlin
  F = mrir_array_accelerated_folding_operator(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for ii = 1:Ncol,

    % transmit matrix, S: Nlin x Ncha
    S = squeeze(transmit(ii, :, :));

    errmtx_inv = zeros(Nlin, Nlin);
    for jj = 1:ceil(Nlin/R),

      % strip of row of folding operator corresponding to this line
      a = F(jj,:);

      % excitation encoding matrix (C: Ncha x Nlin) computes the weighted
      % sum of all "true" pixels in this line for the observed pixel value
      % at each coil. each line is weighted first by the transmit profile of
      % a coil, then the weighted pixels are summed according to the
      % aliasing for this row. "a" is sparse, and the number of nonzero
      % elements is greater than or equal to R.
      C = (diag(a) * S);

      % since in general the pixels contributing to an observation are not
      % unique to the observation, the errors must be summed.
      errmtx_inv = errmtx_inv + (C*C');

    end;

    % desired excitation profile for this column, P: Nlin x 1
    P = reshape(excitation(ii,:), [], 1);
    
    Mu = diag(P);

    %    if ( cond(errmtx_inv) > 1/sqrt(eps) ),
    %      error('the error covariance matrix is singular!');
    %    end;
    errmtx = inv(errmtx_inv);

    % store g-factor for all lines in this image column
    gfactor(ii, :) = sqrt(abs( diag(Mu'*errmtx*Mu) ./ diag([Mu'*inv(diag(diag(errmtx_inv)))*Mu]) ));

  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__general.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
