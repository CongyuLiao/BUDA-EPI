function [epi_dti, epi_dti_raw, epi_dti_img] = mrir_epi_dti(meas)
% MRIR_EPI_DTI
%
%  (wrapper around MRIR_EPI)
  
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2007/dec/13
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  meas.evp.NRepMeas = meas.evp.NSetMeas;
  
  
  [epi_dti, epi_dti_raw, epi_dti_img] = mrir_epi(meas);

  perm = 1:16;
  perm(4) = 7; perm(7) = 4;

  epi_dti     = permute(epi_dti,     perm);
  epi_dti_raw = permute(epi_dti_raw, perm);
  epi_dti_img = permute(epi_dti_img, perm);
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
