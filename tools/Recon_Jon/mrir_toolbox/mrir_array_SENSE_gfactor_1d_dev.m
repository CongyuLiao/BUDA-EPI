function gfactor = mrir_array_Gfactor_1d_dev(raw, noise, R)
%MRIR_ARRAY_GFACTOR_1D_DEV
%
% gfactor = mrir_array_Gfactor_1d_dev(raw, noise, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/06
% $Id: mrir_array_SENSE_gfactor_1d_dev.m,v 1.2 2007/04/06 23:32:21 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % for computations, cast data and noise as doubles to minimize roundoff
  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  sensitivity = mrir_sensitivity_map(double(raw), 100);

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(noise, 'svd');

  sensitivity_whitened = reshape(sensitivity, [], size(noise,3)) * W.';
  sensitivity = reshape(sensitivity_whitened, size(sensitivity));

  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);

  % alias matrix, A: Nlin/R x Nlin
  A = mrir_array_accelerated_aliasop(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for ii = 1:Ncol,

    % sensitivity matrix, S: Nlin x Ncha
    S = squeeze(sensitivity(ii, :, :));

    accumulator = zeros(Nlin, Nlin);
    J = [];
    for jj = 1:Nlin/R,

      % encoding matrix, F: Ncha x Nlin
      F = S.' * diag(A(jj,:))';

      accumulator = accumulator + (F' * F);

      J = cat(1, J, F);

    end;
    noisepower_inv = J'*J;

    %noisepower_inv = accumulator;

    %          E = [];
    noisepower_inv = zeros(Nlin, Nlin);
    for cha=1:Ncha;
      %      E = cat(1, E, A*diag(S(:,cha)));
      E = A*diag(S(:,cha));
      noisepower_inv = noisepower_inv + (E'*E);     
    end;
%    noisepower_inv = E'*E;

    noisepower = inv(noisepower_inv);

    gfactor(ii, :) = sqrt(abs( diag(noisepower) .* diag(noisepower_inv) ));

  end;


  return;

  

%**************************************************************************%
function aliasop = mrir_array_accelerated_aliasop(Nlin, R)

% easy way to calculate aliasing operator: transform full-FOV image indices
% to k-space, remove lines skipped during acceleration, then transform back
% to image space.
  
  fullFOV = mrir_fDFT_freqencode(eye(Nlin));
  reduFOV = fullFOV(1:R:end, :);
  aliasop = mrir_iDFT_freqencode(reduFOV)/Nlin*R;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d_dev.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
