function img_orient = mrir_image_orient_SAG(img)
%MRIR_IMAGE_ORIENT_SAG
%
% img_orient = mrir_image_orient_SAG(img)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/16
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%




  img_orient = flipdim(img, 2);


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
