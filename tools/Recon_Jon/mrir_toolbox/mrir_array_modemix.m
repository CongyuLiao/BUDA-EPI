function varargout = mrir_array_modemix(img, noisecov, varargin)
%MRIR_ARRAY_MODEMIX
%
% mrir_array_modemix(img, noisecov, varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/sep/15
% $Id: mrir_array_modemix.m,v 1.1 2008/10/20 17:48:12 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  if ( nargin >= 3 ),
    M = varargin{1};
  else,

    W = mrir_array_whitening_operator(noisecov, 'svd');
    img_hat = mrir_array_transform(img, W);

    sigcor_hat = mrir_array_stats_matrix(img_hat, 'cor');
    [U_hat, S_hat, V_hat] = svd(sigcor_hat);

    M = W'*U_hat;
  
  end;

  n = size(M, 2);
  if ( nargin >= 4 ),
    n = varargin{2};
  end;

  % truncate
  M = M(:,1:n);
  
  if ( nargout == 1 ),
    varargout{1} = M;
    return;
  end;

  
  %==--------------------------------------------------------------------==%
  %%% apply mode mixing transformation
  
  img_mix = mrir_array_transform(img, M');
  noisecov_mix = M' * noisecov * M;


  if ( nargout > 0 ),

    varargout{1} = img_mix;
    varargout{2} = noisecov_mix;

    varargout{3} = M;
    
    if ( nargin < 3 ),
      varargout{4} = S_hat;
    end;
  
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_modemix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
