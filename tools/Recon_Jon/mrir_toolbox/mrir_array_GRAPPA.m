function img = mrir_array_GRAPPA(meas, varargin)
%

  
% this will be a high-level function that accepts the raw data struct and
% chosen reconstruction parameters, including the kernel size, number of
% reconstruction blocks, and the desired fullFOV image resolution (or
% sampling).

  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/06
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  
  % default parameters for GRAPPA kernel size and number of kernels
  Nx = 3;
  Ny = 3;
  Nz = 3;

  Nk1 = 1;
  Nk2 = 1;

  if ( nargin >= 2 ), Nx  = varargin{1}; end;
  if ( nargin >= 3 ), Ny  = varargin{2}; end;
  if ( nargin >= 4 ), Nz  = varargin{3}; end;
  if ( nargin >= 5 ), Nk1 = varargin{4}; end;
  if ( nargin >= 6 ), Nk2 = varargin{5}; end;


  %==--------------------------------------------------------------------==%

  [dat, acs] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);

  % if ( 2D ),
    
  raw_full = mrir_array_GRAPPA_2d(dat, acs, meas.evp, Nx, Ny, Nz, Nk1, Nk2);

  % end;
  
  
  img = mrir_conventional(raw_full);
  
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
