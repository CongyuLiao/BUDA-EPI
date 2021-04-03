function varargout = epi_GRAPPA_kawin(measfile)
%epi_GRAPPA_kawin
%
% epi_GRAPPA_kawin(measfile)

% NOTE: coordinate Notation is different to the standard one use in pTx need to do [].'
% NOTE: this version is only compatible with combinations that do not use
% any sensitivity weighting!

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/09
% $Id: physio_recon_snr.m,v 1.1 2009/01/09 01:13:22 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  %maxNumCompThreads('automatic');
  
  FLAG__SNR_compute = 1;
  lastwarn('');
  
  [pathstr, filestr] = fileparts(measfile);


  opt.ReturnStruct = 1;
  opt.PhascorCollapseSegments = 0;
  opt.ReadMultipleRepetitions = 0;

  fprintf(1, '\n\n');
  disp(sprintf('reading file "%s"', filestr));

  disp(sprintf(' (repetition %03d) ',   1));
  meas = read_meas_dat__fast(measfile, opt);

  msg = lastwarn;
  if ( strncmp(msg, 'EOF encountered before ACQEND!', 30) ),
    %                                         1 2 3 4 5 6       7 8 9 0 1 2 3 4 5 6
    meas.data =           meas.data(          :,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);
    meas.evp.NRepMeas = mrir_ice_dimensions(meas.data, 'rep');
  end;


  NRep = meas.evp.NRepMeas;
  meas.evp.NRepMeas = 1;

  NCha = meas.evp.NChaMeas;
  NSlc = meas.evp.NSlcMeas;

  evp = meas.evp;
  prot = meas.prot;
  
  % workaround for unusual Siemens convention #1:
  if ( meas.evp.NFirstRefLin == 0 ),
    meas.evp.NFirstRefLin = mrir_ice_dimensions(meas.patrefscan, 'lin') - meas.evp.NRefLin + 1;
  end;


  % workaround for unusual Siemens convention #2:

  % (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
  % last lines are effectively zero-padded; this could throw off SNR
  % calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
  % in same number of k-space lines as there are image lines.)
  if ( ~isempty(evp.NAFLin) && (evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (evp.NLinMeas < evp.NImageLins) ),
    jnotify
    keyboard
    evp.NLinMeas = evp.NImageLins;
  end;

  
  %==--------------------------------------------------------------------==%

  if ( evp.NAFLin > 1 ),
    [dat, acs] = mrir_evi_GRAPPA_prep(meas);         clear meas;
    [epi, G] = mrir_epi_GRAPPA(dat, acs, prot, evp, 3, 4, 1, 1);
  else,
    epi = mrir_epi(meas);                            clear meas;
  end;
  
  epi_rss_reps = mrir_array_combine_sos(epi);


  %==--------------------------------------------------------------------==%

  opt.ReadMultipleRepetitions = 1;
  opt.ExtractRepetition       = 1;


  for rep = 2:NRep,

    disp(sprintf(' (repetition %03d) ', rep));

    meas = read_meas_dat__fast(measfile, opt, rep);

    msg = lastwarn;
    if ( strncmp(msg, 'EOF encountered before ACQEND!', 30) ),
      break;
    end;

    meas.evp  = evp;
    meas.prot = prot;

    if ( evp.NAFLin > 1 ),
      meas.patrefscan = [];
      meas.patrefscan_phascor = [];
      
      dat = mrir_evi_GRAPPA_prep(meas);                  clear meas;
 
% Kawin 3/24/2009: change so that it will read in multiple measurements.....
      epi(:,:,:,rep) = mrir_epi_GRAPPA(dat, G, prot, evp, 3, 4, 1, 1);
%      epi = mrir_epi_GRAPPA(dat, G, prot, evp, 3, 4, 1, 1);
    else,    
      epi(:,:,:,rep) = mrir_epi(meas);                              clear meas;
    end;
    
  end;


  %==--------------------------------------------------------------------==%

  if ( nargout > 0 ),
     % Kawin 3/24/2009: change so that it will read in multiple measurements.....
     % varargout{1} = sum(epi,4)/NRep;
     epi = conj( -i*epi);%this make it the same notation as what I have in ReadMeasDotDat
     varargout{1} = epi; 
     eval(['save ' measfile(1:end-4) '.mat epi' ]);
  end;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/physio_toolbox/physio_recon_snr.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
