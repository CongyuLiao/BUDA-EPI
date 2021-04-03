function img_smoothed = mrir_image_smooth_3d(img, fwhm_voxel)
%MRIR_IMAGE_SMOOTH_3D
%
% img_smoothed = mrir_image_smooth_3d(img, fwhm_voxel)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/07
% $Id: mrir_image_smooth_3d.m,v 1.1 2008/11/06 06:07:45 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  sigma_voxel = fwhm_voxel / ( sqrt(2*log(2))*2 )

  dims = size(img);

  [x,y,z] = meshgrid(linspace(-dims(1)/2, +dims(1)/2, dims(1)), linspace(-dims(2)/2, +dims(2)/2, dims(2)), linspace(-dims(9)/2, +dims(9)/2, dims(9)));

  h = exp(-(x.*x + y.*y + z.*z)/(2*sigma_voxel^2));
  h = h / sum(h(:));

  H = mrir_fDFT(mrir_fDFT(mrir_fDFT(  h, 1, dims(1)), 2, dims(2)), 3, dims(9));


  img_smoothed = [];
  for cha = 1:mrir_ice_dimensions(img, 'cha'),

    img_cha = squeeze(img(:,:,cha, 1,1,1,1,1, :));

    D = mrir_fDFT(mrir_fDFT(mrir_fDFT(img_cha, 1), 2), 3);
    S = D .* H;

    img_smoothed(:,:,cha,:) = mrir_iDFT(mrir_iDFT(mrir_iDFT(S, 3), 2), 1);

  end;

  img_smoothed = reshape(img_smoothed, size(img));


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_smooth_3d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
