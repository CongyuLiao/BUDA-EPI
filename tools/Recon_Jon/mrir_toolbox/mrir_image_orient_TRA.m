function img_orient = mrir_image_orient_TRA(img)
%MRIR_IMAGE_ORIENT_TRA
%
% img_orient = mrir_image_orient_TRA(img)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/16
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims = 1:ndims(img);
  dims_perm = [2 1 dims(3:end)];

  img_orient = flipdim(permute(img, dims_perm),1);


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
