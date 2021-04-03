function img_dilate = mrir_volume_dilate(img, radius)

% img_dilate = mrir_volume_dilate(img, radius)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  [xx, yy, zz] = meshgrid(-radius:+radius);
  sphere = sqrt(xx.^2 + yy.^2 + zz.^2) <= radius;


  img_dilate = [];
  for cha = 1:mrir_ice_dimensions(img, 'cha'),

    img_chan = squeeze(img(:,:,cha, 1,1,1,1,1, :));

    if ( isreal(img_chan) ),
      img_dilate(:,:,cha,:) = imdilate(img_chan, strel('arbitrary', sphere));
    else,
      
      img_dilate_real = imdilate(real(img_chan), strel('arbitrary', sphere));
      img_dilate_imag = imdilate(imag(img_chan), strel('arbitrary', sphere));
      
      img_dilate(:,:,cha,:) = complex(img_dilate_real, img_dilate_imag);
      
    end;

   
  end;

  img_dilate = reshape(img_dilate, size(img));
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

