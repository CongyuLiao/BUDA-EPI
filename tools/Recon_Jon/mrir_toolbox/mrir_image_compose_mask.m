function img_composite = mrir_image_compose_mask(img1, img0, mask)
%MRIR_IMAGE_COMPOSE_MASK  compose two images with a binary mask
%
% img_composite = mrir_image_compose_mask(img1, img0, mask)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/oct/20
% $Id: mrir_image_compose_mask.m,v 1.1 2008/10/21 03:44:35 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims_imgs = size(img1);
  dims_mask = size(mask);

  dims_mask = [dims_mask, ones(1, ndims(img1) - ndims(mask))];


  mask_rep = repmat(mask, dims_imgs./dims_mask);


  img_composite = [img1 .* (mask_rep)] + [img0 .* (~mask_rep)];


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_compose_mask.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
