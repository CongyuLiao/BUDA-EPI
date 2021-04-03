function img_smooth = mrir_image_antialias(img)
%MRIR_IMAGE_ANTIALIAS  apply anti-aliasing smoothing akin to SYNGO's "IPT"
  
%
% img_smooth = mrir_image_antialias(img)
    
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/21
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  for ind = 1:size(img,3),
  
  img_smooth(:,:,ind) = imresize(img(:,:,ind), 2, 'nearest');
  
  
  end;
  
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
