function condition = mrir_array_SENSE_condition_number(sensitivity, covmtx, ...
						  R1, R2)
%MRIR_ARRAY_SENSE_CONDITION_NUMBER  condition number of SENSE reconstruction
%
% condition = mrir_array_SENSE_condition_number(sensitivity, covmtx, R1, R2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/11
% $Id: mrir_array_condition_number.m,v 1.3 2007/04/06 22:21:31 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtxinv = inv(covmtx);

  aliasmap1 = mrir_array_accelerated_aliasmap(Nlin, R1);
  aliasmap2 = mrir_array_accelerated_aliasmap(Ncol, R2);

  if ( prod(size(aliasmap1)) ~= Nlin )
    warning('R1 not compatible with Nlin');
  end;
 
  if ( prod(size(aliasmap2)) ~= Ncol )
    warning('R2 not compatible with Ncol');
  end;
  
  % preallocate
  condition  = zeros(Ncol, Nlin);

  for jj = 1:floor(Nlin/R1),

    % indices of the aliased pixels for this position from lookup table
    aliased_ind1 = aliasmap1(jj, 1:R1);

    for ii = 1:floor(Ncol/R2),

      % indices of the aliased pixels for this position from lookup table
      aliased_ind2 = aliasmap2(ii, 1:R2);

      % R2 x R1 x Nchan
      b = sensitivity(aliased_ind2, aliased_ind1, :);

      % encoding matrix: Nchan x R1*R2, or Nchan x Nalias
      E = permute(reshape(b, [R1*R2, Ncha]), [2, 1]);

      % information matrix
      nu = E' * covmtxinv * E;

      condition(aliased_ind2,aliased_ind1) = cond(nu);

    end;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/COIL_ARRAYS/MATLAB/mrir_array_condition_number.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
