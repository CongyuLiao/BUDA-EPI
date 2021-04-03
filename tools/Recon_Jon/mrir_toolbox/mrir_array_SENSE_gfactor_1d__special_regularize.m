function gfactor = mrir_array_SENSE_gfactor_1d__special(sensitivity, covmtx, R, varargin)
%MRIR_ARRAY_SENSE_GFACTOR_1D__SPECIAL_REGULARIZE
%
% gfactor = mrir_array_SENSE_gfactor_1d__special(sensitivity, covmtx, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id: mrir_array_Gfactor_1d__special.m,v 1.1 2007/05/05 03:35:28 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  f = 0.0000005;

  if ( nargin >= 4 ),
    f = varargin{1};
  end;
  
  
  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtxinv = inv(covmtx);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin);
  artifact = zeros(Ncol, Nlin);

  for jj = 1:floor(Nlin/R),

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);

    for ii = 1:Ncol,
      % encoding matrix, E: Ncha x R
      E = permute(squeeze( sensitivity(ii, aliased_pixels_ind, :) ), [2 1]);

      errcovmtx_inv = E' * covmtxinv * E;

      if ( cond(errcovmtx_inv) > 1/sqrt(eps) ),
        error('the error covariance matrix is singular!');
      end;

      e = (max(eig(errcovmtx_inv)));

      % kellman2001SENSE
      % sodickson2000tailored
      
      errcovmtx = inv( errcovmtx_inv + f*e*eye(size(errcovmtx_inv)) );

      U = errcovmtx * E' * covmtxinv;
      rho = U * E;
      A = diag(1./diag(rho));

      nu = A * rho * errcovmtx * A;
      
      
      gfactor(ii, aliased_pixels_ind) = sqrt(abs( diag(nu) .* diag(errcovmtx_inv) ));

      %rho = A * rho;
      %artifact(ii, aliased_pixels_ind) = [rho(1,2), rho(2,1)];
      

    end;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_Gfactor_1d__special.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: