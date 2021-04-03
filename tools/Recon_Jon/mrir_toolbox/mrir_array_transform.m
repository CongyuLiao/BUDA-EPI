function S_hat = mrir_array_transform(S, W)
%MRIR_ARRAY_TRANSFORM
%
% S_hat = mrir_array_transform(S, W)

% (formerly "mrir_array_whitening_apply.m")
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jul/19
% $Id: mrir_array_transform.m,v 1.1 2008/04/26 01:16:48 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  
  
  dims = size(S);
  dims(end+1:16) = 1;

  Ncha = dims(3);

  dims(3) = size(W, 1);
  
  Nrem = prod(dims([1:2,4:end]));
  
  % permute and reshape n-dimension array to Ncha x Nrem matrix
  s = reshape(permute(S, [3, 1, 2, 4:16]), Ncha, []);
  
  if ( (size(s,1) ~= Ncha) || (size(s,2) ~= Nrem) ),
    error('incorrect dimensions found in resized matrix');
  end;
  
  s_hat = W * s;
  
  S_hat = ipermute(reshape(s_hat, dims([3, 1, 2, 4:16])), [3, 1, 2, 4:16]);
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_transform.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
