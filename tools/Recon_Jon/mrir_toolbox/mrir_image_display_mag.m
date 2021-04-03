function varargout = mrir_image_display_mag(data_complex)
%MRIR_IMAGE_DISPLAY_MAG
%
% [ph, ax, fig] = mrir_image_display_mag(data_complex)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/07
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims = size(data_complex);

  fig = figure;
  ax = axes;

  ph = imagesc(abs(data_complex).');
  axis image;
  colormap(gray);

  set(ax, 'DataAspectRatio', [1 1 1]);

  if ( nargout > 0 ),
    varargout{1} = ph;
    varargout{2} = ax;
    varargout{3} = fig;
  end;


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
