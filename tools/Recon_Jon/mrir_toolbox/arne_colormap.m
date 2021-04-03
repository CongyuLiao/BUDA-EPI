function varargout = arne_colormap(varargin)


  cmap = [0            0      0.88547
          0            0            1
          0      0.71922            1
          0      0.88547            1
          0            1            1
          0.71922            1      0.88547
          0.88547            1      0.71922
          1            1            0
          1      0.88547            0
          1      0.71922            0];

  
  if ( nargout >= 1 ),
    varargout{1} = cmap;
  else,
    colormap(cmap);
    cb = colorbar;
    set(cb, 'YTick', 0.1:0.1:1.0)
  end;
  
  
  return;