function meas = mrir_measdat_check(meas)
%MRIR_MEASDAT_CHECK
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/sep/20
% $Id: mrir_measdat_check.m,v 1.2 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [meas, FLAG__fix_7T_channels] = mrir_hack__fix_7T_channels(meas);

  if ( FLAG__fix_7T_channels ),
    disp(sprintf('<i> [%s]: repaired meas.dat file---Bay 5 RCCS channel assignment issue with 32-channel coil FIXED', mfilename));
  end;

  msg = lastwarn;
  if ( strcmp(meas.STATUS, 'ABORTED') ),

    disp(sprintf('==> [%s]: repairing meas.dat file---discarding final repetition of aborted scan...', mfilename));

    %                                         1 2 3 4 5 6       7 8 9 0 1 2 3 4 5 6
    meas.data =           meas.data(          :,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);
    meas.evp.NRepMeas = mrir_ice_dimensions(meas.data, 'rep');

    disp(sprintf('<i> [%s]: truncation complete.', mfilename));
    meas.STATUS = 'TRUNCATED';

  end;


  NRep_read = meas.evp.NRepMeas;
  NRep_meas = NRep_read;
  meas.evp.NRepMeas = 1;

  % workaround for unusual Siemens convention #1:
  if ( meas.evp.NFirstRefLin == 0 & isfield(meas, 'patrefscan') ),
    meas.evp.NFirstRefLin = mrir_ice_dimensions(meas.patrefscan, 'lin') - meas.evp.NRefLin + 1;
  end;

  if ( meas.evp.NFirstRefPar == 0 & isfield(meas, 'patrefscan') ),
    meas.evp.NFirstRefPar = mrir_ice_dimensions(meas.patrefscan, 'par') - meas.evp.NRefPar + 1;
  end;

  
  % workaround for unusual Siemens convention #2:

  % (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
  % last lines are effectively zero-padded; this could throw off SNR
  % calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
  % in same number of k-space lines as there are image lines.)
  if ( ~isempty(meas.evp.NAFLin) && (meas.evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (meas.evp.NLinMeas < meas.evp.NImageLins) ),
    keyboard
    meas.evp.NLinMeas = meas.evp.NImageLins;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_measdat_check.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
