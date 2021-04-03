function [evi, evi_raw, evi_img] = mrir_evi_patrefscan(meas)
%MRIR_EVI_PATREFSCAN
%
% [evi, evi_raw, evi_img] = mrir_evi_patrefscan(meas)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  meas_patrefscan = meas;

  meas_patrefscan.data           = meas.patrefscan;
  meas_patrefscan.data_phascor1d = meas.patrefscan_phascor;

  meas_patrefscan = rmfield(meas_patrefscan, {'patrefscan', 'patrefscan_phascor'});

  % assume no partial fourier on the patrefscan lines
  meas_patrefscan.prot.ucPhasePartialFourier = 1;

  % assume only one repetition for EVI patrefscan
  meas_patrefscan.evp.NRepMeas = 1;

  [evi, evi_raw, evi_img] = mrir_evi(meas_patrefscan);


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
