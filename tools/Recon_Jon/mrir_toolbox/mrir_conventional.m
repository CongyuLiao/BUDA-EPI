function img_uncomb = mrir_conventional(raw, varargin)
%MRIR_CONVENTIONAL  reconstructs conventional (cartesian) acquisitions
%
% img_uncomb = mrir_conventional(raw, prot);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  prot = [];
  if ( nargin >= 2 ), prot = varargin{1}; end;
  if ( isempty(prot) ), prot = read_meas_prot__struct; end;
 
  
  %==--------------------------------------------------------------------==%

  % average k-space lines across averages
  raw_avg = mean(raw, 16);
  
  
  Npar = size(raw_avg, 9);
  
  if ( Npar == 1 ),
    img_uncomb = mrir_conventional_2d(raw_avg, prot);
  else,
    img_uncomb = mrir_conventional_3d(raw_avg, prot);
  end;

  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
