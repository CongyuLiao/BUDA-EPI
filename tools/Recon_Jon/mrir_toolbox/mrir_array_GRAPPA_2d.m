function raw = mrir_array_GRAPPA_2d(dat, acs, evp, varargin)
%MRIR_ARRAY_GRAPPA_2D
%
% raw = mrir_array_GRAPPA_2d(dat, acs, evp)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/14
% $Id: mrir_array_GRAPPA_2d.m,v 1.2 2008/11/22 20:10:44 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
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

  G = mrir_array_GRAPPA_2d_kernel(acs, ...
                                  evp.NAFLin, evp.NAFPar, ...
                                  Nx, Ny, Nz, Nk1, Nk2);

  for rep = 1:evp.NRepMeas,
    for slc = 1:evp.NSlcMeas,
      
      %   1 2 3 4 5 6    7 8 9    0 1 2 3 4 5 6
      raw(:,:,:,:,:,:, rep,:,:, slc,:,:,:,:,:,:) ...
	  = mrir_array_GRAPPA_2d_recon(...
	      dat(:,:,:,:,:,:, rep,:,:, slc,:,:,:,:,:,:), G, ...
	      evp.NAFLin, evp.NAFPar, ...
	      evp.NLinMeas, evp.NParMeas, ...
	      evp.NFirstLin, evp.NFirstPar);
    end;
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
