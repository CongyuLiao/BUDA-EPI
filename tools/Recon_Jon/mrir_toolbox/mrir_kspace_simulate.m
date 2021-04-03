function k_raw = mrir_kspace_simulate(img)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/11
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(img, 1);
  Nlin = size(img, 2);
    
  img_uncrop = padarray(img, [Ncol/2, 0], 'replicate', 'both');

  [ii,jj] = meshgrid(1:Nlin, 1:Ncol*2);

  
  Ncycles = 5;
  
  x_img = cos(2*pi*ii/(Ncol*2)*Ncycles) .* img_uncrop;
  y_img = sin(2*pi*ii/(Ncol*2)*Ncycles) .* img_uncrop;

  z_img = x_img + i*y_img;
  
  z_img = z_img + complex(randn(size(z_img)), randn(size(z_img)))*0.05;
  
  
  
  k_raw = mrir_fDFT_freqencode(mrir_fDFT_phasencode(z_img));
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
