function img_combine_cov = mrir_array_combine_cov(img_uncombined, Rn)
%MRIR_ARRAY_COMBINE_COV
%
% img_combine_cov = mrir_array_combine_cov(img_uncombined, Rn)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/26
% $Id: mrir_array_combine_cov.m,v 1.1 2008/01/27 00:27:29 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
 
  Ncol = mrir_ice_dimensions(img_uncombined, 'col');
  Nlin = mrir_ice_dimensions(img_uncombined, 'lin');
  Ncha = mrir_ice_dimensions(img_uncombined, 'cha');
  
  % arbitrary constant (normalizes image to approx same level as RSS combo)
  %alpha = sqrt(mean(diag(Rn)));
  alpha = 1; % Kawin Setsompop: to make things easier for grappa analytic noise calc etc
  
  dims = size(img_uncombined);
  dims(end+1:16) = 1;
  
  
  img_permute = permute(img_uncombined, [3, 1, 2, 4:ndims(img_uncombined)]);
  img_reshape = reshape(img_permute, Ncha, []);
  
  img_reshape_cov = zeros(1, size(img_reshape, 2));

  iRn = inv(Rn);
  
  for ind = 1:length(img_reshape_cov),
  
    img_reshape_cov(ind) = abs ( alpha * sqrt(abs(img_reshape(:,ind)' * iRn * img_reshape(:,ind))) );

  end;
    
  img_combine_cov = reshape(img_reshape_cov, [Ncol, Nlin, 1, dims(4:end)]);
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_cov.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
