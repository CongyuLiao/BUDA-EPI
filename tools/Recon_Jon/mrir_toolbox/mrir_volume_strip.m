function img_strip = mrir_volume_strip(img, radius)

% img_strip = mrir_volume_strip(img, radius)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  [xx, yy, zz] = meshgrid(-radius:+radius);
  sphere = sqrt(xx.^2 + yy.^2 + zz.^2) <= radius;


  img_strip = [];
  for cha = 1:mrir_ice_dimensions(img, 'cha'),

    img_chan = squeeze(img(:,:,cha, 1,1,1,1,1, :));

    if ( isreal(img_chan) ),
      img_strip(:,:,cha,:) = imerode(img_chan, strel('arbitrary', sphere));
    else,
      
      img_strip_real = imerode(real(img_chan), strel('arbitrary', sphere));
      img_strip_imag = imerode(imag(img_chan), strel('arbitrary', sphere));
            
      img_strip(:,:,cha,:) = complex(img_strip_real, img_strip_imag);
      
    end;

  end;

  img_strip = reshape(img_strip, size(img));

  
  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

