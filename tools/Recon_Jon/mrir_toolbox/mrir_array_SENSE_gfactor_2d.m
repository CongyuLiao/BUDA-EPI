function gfactor = mrir_array_Gfactor_2d(raw, noise, R1, R2)
%MRIR_ARRAY_GFACTOR_2D
%
% gfactor = mrir_array_Gfactor_2d(raw, noise, R1, R2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/04
% $Id: mrir_array_SENSE_gfactor_2d.m,v 1.4 2007/06/20 01:45:20 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % for computations, cast data and noise as doubles to minimize roundoff
  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  sensitivity = mrir_sensitivity_map(double(raw));

  covmtxinv = inv(covmtx);

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(noise, 'svd');

  sensitivity_whitened = reshape(sensitivity, [], size(noise,3)) * W.';
  sensitivity = reshape(sensitivity_whitened, size(sensitivity));


  Ncol = size(sensitivity, 1);
  Nlin = size(sensitivity, 2);
  Ncha = size(sensitivity, 3);

  alias_map1 = mrir_array_accelerated_aliasmap(Nlin, R1);
  alias_map2 = mrir_array_accelerated_aliasmap(Ncol, R2);

  if ( prod(size(alias_map1)) ~= Nlin )
    warning('R1 not compatible with Nlin');
  end;
 
  if ( prod(size(alias_map2)) ~= Ncol )
    warning('R2 not compatible with Ncol');
  end;
  
  % preallocate
  gfactor  = zeros(Ncol, Nlin);


  for jj = 1:floor(Nlin/R1),

    % indices of the aliased pixels for this position from lookup table
    aliased_ind1 = alias_map1(jj, 1:R1);

    for ii = 1:floor(Ncol/R2),

      % indices of the aliased pixels for this position from lookup table
      aliased_ind2 = alias_map2(ii, 1:R2);

      % R2 x R1 x Nchan
      b = sensitivity(aliased_ind2, aliased_ind1, :);

      % encoding matrix: Nchan x R1*R2, or Nchan x Nalias
      E = permute(reshape(b, [R1*R2, Ncha]), [2, 1]);

      % information matrix (whitened)
      nu = E' * E;
      eta = inv(nu);

      gfactor(aliased_ind2,aliased_ind1) = reshape(sqrt(abs(diag(nu).*diag(eta))), [R2,R1]);

    end;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
