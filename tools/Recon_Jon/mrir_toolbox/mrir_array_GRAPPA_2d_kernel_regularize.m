function varargout = mrir_array_GRAPPA_2d_kernel_regularize(k_ACS, R1, R2, varargin)
%MRIR_ARRAY_GRAPPA_2D_KERNEL_REGULARIZE
%
% G = mrir_array_GRAPPA_2d_kernel_regularize(k_ACS, R1, R2)
% G = mrir_array_GRAPPA_2d_kernel_regularize(k_ACS, R1, R2, Nx, Ny, Nz)
% G = mrir_array_GRAPPA_2d_kernel_regularize(k_ACS, R1, R2, Nx, Ny, Nz, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2008/may/01
% $Id: mrir_array_GRAPPA_2d_kernel.m,v 1.3 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  Nx = 3;
  Ny = 3;
  Nz = 3;
  Nk1 = 1;
  Nk2 = 1;

  if ( nargin >= 4 ),
    Nx  = varargin{4-3};
  end;

  if ( nargin >= 5 ),
    Ny  = varargin{5-3};
  end;

  if ( nargin >= 6 ),
    Nz  = varargin{6-3};
  end;

  if ( nargin >= 7 ),
    Nk1  = varargin{7-3};
  end;

  if ( nargin >= 8 ),
    Nk2  = varargin{8-3};
  end;


  %==--------------------------------------------------------------------==%


  global REGULARIZE; if ( isempty(REGULARIZE) ), REGULARIZE = 1; end;


  varargout = mrir_array_GRAPPA_2d_kernel(k_ACS, R1, R2, Nx, Ny, Nz, Nk1, Nk2);




  clear REGULARIZE


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_kernel.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
