function varargout = mrir_sensitivity_map__estimate_update_leo(img_full, img_redu, sens_initial, covmtx)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/23
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%





  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:





% On Mon, 30 Apr 2007, Grady, Leo (SCR US) wrote:
%
% > Hi Jon,
% >
% > Thanks for your help earlier.
% >
% >
% > Just to follow up on the mathematics here, consider a problem:
% > E \tilde{x} = y
% > where you solve for \tilde{x}, given an inaccurate E but an accurate y.
% >
% > Now, assume that you obtain, somehow, an accurate value for x and you
% > want to find a minimal update U, in order to produce E^* = E+U, such
% > that
% > E^* x = y
% >
% > We set up the problem like this:
% > min_U trace(U^T U)
% > s.t. (E+U) x = y
% >
% > After some variational calculus and algebra, you find that
% > U = - \frac{E x x^T - y x^T} {x^T x}
% >
% >
% > Let me know what you think.
% >
% >
% > Leo
% >
