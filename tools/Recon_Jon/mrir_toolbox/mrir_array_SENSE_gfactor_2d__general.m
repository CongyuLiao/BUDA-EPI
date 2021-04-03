function [gfactor, condition] = mrir_array_SENSE_gfactor_2d__general(sensitivity, covmtx, R1, R2)
%MRIR_ARRAY_SENSE_GFACTOR_2D__GENERAL
%
% 
%  DOES NOT WORK!
%
% gfactor = mrir_array_Gfactor_1d(raw, noise, R1, R2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/04
% $Id: mrir_array_SENSE_gfactor_2d__general.m,v 1.1 2007/04/24 19:48:42 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  disp('not yet!');
  return;
  
  % for computations, cast data and noise as doubles to minimize roundoff
%  sensitivity = mrir_sensitivity_map(double(raw), 100);

  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(covmtx, 'svd');

  sensitivity_hat = mrir_array_transform(sensitivity_hat, W);

  % alias matrix, A: ceil(Nlin/R) x Nlin
  A1 = mrir_array_accelerated_aliasop(Nlin, R1);
  A2 = mrir_array_accelerated_aliasop(Ncol, R2);

  % preallocate
  gfactor  = zeros(Ncol, Nlin);

  for jj = 1:Nlin/R1,

    sensitivity_folded = squeeze(sensitivity(:, jj, :))

  end;


  for ii = 1:Ncol,

    % sensitivity matrix, S: Nlin x Ncha
    S = squeeze(sensitivity(ii, :, :));

    errcovmtx_inv = zeros(Nlin, Nlin);
    for channel = 1:Ncha;

      % per-channel encoding matrix, E: ceil(Nlin/R) x Nlin

      % scale the columns of A by the column of S for this channel
      E = A * diag( S(:, channel) );

      % accumulate the inner product across channels (alternatively, could
      % just vertically stack E into new matrix F: (Ncha*Nlin/R) x Nlin so
      % that errcovmtx_inv = F'*F.)
      errcovmtx_inv = errcovmtx_inv + (E'*E);
    end;

    if ( cond(errcovmtx_inv) > 1/sqrt(eps) ),
      error('noise power matrix is singular!');
    end;
    errcovmtx = inv(errcovmtx_inv);

    % store g-factor for all lines in this image column
    gfactor(ii, :) = sqrt(abs( diag(errcovmtx) .* diag(errcovmtx_inv) ));

  end;


  return;



%**************************************************************************%
function aliasop = mrir_array_accelerated_aliasop(N, R)

% easy way to calculate aliasing operator: transform full-FOV image indices
% to k-space, remove lines skipped during acceleration, then transform back
% to image space.

  fullFOV = mrir_fDFT_freqencode(eye(N));
  reduFOV = fullFOV(1:R:end, :);
  aliasop = mrir_iDFT_freqencode(reduFOV)/N*R;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/COIL_ARRAYS/MATLAB/mrir_array_SENSE_gfactor_2d__general.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
