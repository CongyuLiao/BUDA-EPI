function gfactor = mrir_array_Gfactor_1d(raw, noise, R)
%MRIR_ARRAY_GFACTOR_1D
%
% gfactor = mrir_array_Gfactor_1d(raw, noise, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/11
% $Id: mrir_array_Gfactor_1d.m,v 1.3 2007/04/06 22:21:31 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Nlin = size(raw, 2); % number of phase encoded lines

  % test if the number of phase encoded lines are an integer multiple of the
  % acceleration R
  if ( mod(Nlin, R) == 0 ),

    % since Nlin divisible by R, we can exploit a faster method for this special case
    gfactor = mrir_array_Gfactor_1d__special(raw, noise, R);
  else,

    % Nlin not divisible by R: must use slower algorithm
    disp('number of PE lines not integer multiple of R! using slower, general method...');
    gfactor = mrir_array_Gfactor_1d__general(raw, noise, R);
  end;


  return;



%**************************************************************************%
function gfactor = mrir_array_Gfactor_1d__special(raw, noise, R)

  % for computations, cast data and noise as doubles to minimize roundoff
  sensitivity = mrir_sensitivity_map(double(raw), 100);

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  covmtxinv = inv(covmtx);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for jj = 1:Nlin/R,

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);

    for ii = 1:Ncol,
      % encoding matrix, E: Ncha x R
      E = permute(squeeze( sensitivity(ii, aliased_pixels_ind, :) ), [2 1]);

      noisepower_inv = E' * covmtxinv * E;

      if ( cond(noisepower_inv) > 1/sqrt(eps) ),
	error('noise power matrix is singular!');
      end;
      noisepower = inv(noisepower_inv);

      gfactor(ii, aliased_pixels_ind) = sqrt(abs( diag(noisepower) .* diag(noisepower_inv) ));

    end;
  end;


  return;



%**************************************************************************%
function gfactor = mrir_array_Gfactor_1d__general(raw, noise, R)

  % for computations, cast data and noise as doubles to minimize roundoff
  sensitivity = mrir_sensitivity_map(double(raw), 100);

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(noise, 'svd');

  sensitivity_whitened = reshape(sensitivity, [], size(noise,3)) * W.';
  sensitivity = reshape(sensitivity_whitened, size(sensitivity));

  % alias matrix, A: ceil(Nlin/R) x Nlin
  A = mrir_array_accelerated_aliasop(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for ii = 1:Ncol,

    % sensitivity matrix, S: Nlin x Ncha
    S = squeeze(sensitivity(ii, :, :));

    % encoding matrix, E: Nlin x Nlin
    E = (A * S).';
    
    E = E * A;
    
    
    noisepower_inv = E'*E;

    if ( cond(noisepower_inv) > 1/sqrt(eps) ),
      error('noise power matrix is singular!');
    end;
    noisepower = inv(noisepower_inv);

    % store g-factor for all lines in this image column
    gfactor(ii, :) = sqrt(abs( diag(noisepower) .* diag(noisepower_inv) ));

  end;


  return;



%**************************************************************************%
function [aliasmap, varargout] = mrir_array_accelerated_aliasmap(Nlin, R, ind_line)

  % "aliasmap_full2redu" is a mapping from full-FOV pixel indices to
  % reduced-FOV pixel indices caused by subsampling during acceleration; and
  % "aliasmap" is the inverse mapping.

  % aliasmap_full2redu:      1 x Nlin  (each FULL maps to one REDU)
  % aliasmap          : Nlin/R x R     (each REDU maps to R FULL)

  % shift map by half a period if R is even
  if ( mod(R, 2) == 0 ),
    shift = floor(floor(Nlin/R) / 2);
  else,
    shift = 0;
  end;

  % the mod operator provides the aliasing mapping
  aliasmap_full2redu = mod( [1:Nlin] + shift, floor(Nlin/R));
  aliasmap_full2redu(find(aliasmap_full2redu==0)) = floor(Nlin/R);

  % invert table to map reduced FOV indices to full FOV indices

  % (since Nlin might not be evenly divisible by R, only find the first R matches)
  for ind = 1:(floor(Nlin/R)),
    aliasmap(ind, 1:R) = find(aliasmap_full2redu == ind, R);
  end;

  % calculate and return aliasop if requested
  if ( nargout >= 2 ),
    aliasop = zeros(R, Nlin);
    aliasop(sub2ind([R, Nlin], 1:R, aliasmap(ind_line,:))) = 1;
    varargout{1} = aliasop;
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

%  aliasop = [zeros( (Nlin - ceil(Nlin/R))/2, Nlin); 
%	     aliasop; 
%	     zeros( (Nlin - ceil(Nlin/R))/2, Nlin)];

  
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/COIL_ARRAYS/MATLAB/mrir_array_Gfactor_1d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
