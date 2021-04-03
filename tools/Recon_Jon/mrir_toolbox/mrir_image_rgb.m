function img_rgb = mrir_image_rgb(img, varargin)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/30
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  irange = [];

  N = 64;
  cmap = jet(64);


  if ( nargin >= 2 ),
    irange = varargin{1};
  end;

  if ( nargin >= 3 ),
    cmap = varargin(2);
    N = size(cmap, 1);
  end;


  %==--------------------------------------------------------------------==%

  if ( isempty(irange) ),
    %                                     m        g   i      r
    img_rgb = ind2rgb(gray2ind(mat2gray(img        ), N), cmap);
  else,
    img_rgb = ind2rgb(gray2ind(mat2gray(img, irange), N), cmap);
  end;

  
  img_r = img_rgb(:,:,1);
  img_g = img_rgb(:,:,2);
  img_b = img_rgb(:,:,3);
  
  img_r(~isfinite(img)) = 0;
  img_g(~isfinite(img)) = 0;
  img_b(~isfinite(img)) = 0;
  
  img_rgb = cat(3, img_r, img_g, img_b);
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
