function varargout = mrir_image_display_aspect(Mvxl2lph, dims, varargin)
%MRIR_IMAGE_DISPLAY_ASPECT
%
% aspect = mrir_image_display_aspect(Mvxl2lph, dims)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/07
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  aspect = [norm(Mvxl2lph(:,1)), norm(Mvxl2lph(:,2)), norm(Mvxl2lph(:,3))];

  if ( nargin >= 3 ),
    handle = varargin{1};

    if ( strcmp(get(handle, 'Type'), 'figure') ),
      AX = min(findall(handle, 'Type', 'axes'));
    elseif ( strcmp(get(handle, 'Type'), 'axes') ),
      AX = handle;
    else
      error('graphics handle must be figure or axes');
    end;

  else,
    AX = gca;
  end;


  set(AX, 'DataAspectRatio', [aspect(dims(1)), aspect(dims(2)), 1]);

  if ( nargout > 0 ),
    varargout{1} = aspect;
  end;


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
