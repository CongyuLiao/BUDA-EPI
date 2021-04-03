function evi_img = mrir_evi_GRAPPA_ice(meas, varargin)
%MRIR_EVI_GRAPPA_ICE
%
% evi_img = mrir_evi_GRAPPA_ice(meas)
% evi_img = mrir_evi_GRAPPA_ice(meas, Nx, Ny, Nz)
% evi_img = mrir_evi_GRAPPA_ice(meas, Nx, Ny, Nz, Nk1, Nk2)
%
%
% see also MRIR_EVI_GRAPPA, MRIR_EVI_GRAPPA_PREP, MRIR_EVI_GRAPPA_PREP_ICE.


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/nov/14
% $Id: mrir_evi_GRAPPA_ice.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
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

  [dat, acs] = mrir_evi_GRAPPA_prep(meas);

  %                           1 2 3 4 5 6 7 8                   9 0 1 2 3 4 5 6
%  acs = acs(end*0.25+1:end*0.75,:,:,:,:,:,:,:,end*0.25+1:end*0.75,:,:,:,:,:,:,:);
  
  evi_img = mrir_evi_GRAPPA(dat, acs, meas.prot, meas.evp, Nx, Ny, Nz, Nk1, Nk2);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_GRAPPA_ice.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
