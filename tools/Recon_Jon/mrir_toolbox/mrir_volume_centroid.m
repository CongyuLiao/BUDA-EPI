function centroid = mrir_volume_centroid(img)
  
% DOES NOT WORK

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % how can we pass the arguments to ndgrid?
  grid = ndgrid(deal(size(img)));
  
  
  centroid_ind = sum(grid(:).*img(:)) ./ sum(img(:));
  
  
  centroid = ind2sub(size(img), centroid_ind);
  
  
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
