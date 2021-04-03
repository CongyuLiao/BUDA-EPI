function [epi, epi_ksp, epi_img] = mrir_epi_debug(meas)
%MRIR_EPI
%
% epi_img = mrir_epi_debug(meas_struct)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/16
% $Id: mrir_epi.m,v 1.5 2007/11/17 00:47:38 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  global FLAG__FIRST_CALL;
  FLAG__FIRST_CALL = 1;
  
  for rep = 1:meas.evp.NRepMeas,
    for slc = 1:meas.evp.NSlcMeas,
      
      % siemens: (1) RO regrid, (2) RO apodize, (3) ROFT, (4) PE correct, (5) PE apod, (6) PEFT

      % extract data
      epi_raw = meas.data(:,:,:,1,1,1,rep,:,1,slc);

      % extract navigators for 1D phase correction
      phascor1d = meas.data_phascor1d(:,:,:,1,1,1,rep,:,1,slc);
      
      
      %==----------------------------------------------------------------==%
      %%% {1}: regridding to correct for ramp sampling
 
      % skip the regridding step if neither trapezoid or sinusoid
      if ( meas.prot.alRegridMode > 1 ),
	if ( FLAG__FIRST_CALL ),
	  prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(epi_raw, 1));
	end;
      
	% regrid *before* phase correction, and deapodize
	epi_raw_grid = mrir_regrid_trapezoid_alt(epi_raw, meas.prot, prot_trapezoid, 1);

	phascor1d_grid = mrir_regrid_trapezoid_alt(phascor1d, meas.prot, prot_trapezoid, 1);
	
      else,	
	% skip
	epi_raw_grid = epi_raw;
      end;


      %==----------------------------------------------------------------==%
      %%% {2}: 1-D apodization filter to reduce ringing due to truncation

      epi_raw_apod = mrir_filter_raw_apodize_1d(epi_raw_grid, 1, 0.15);


      %==----------------------------------------------------------------==%
      %%% {3}: iDFT to calculate hybrid (both k space and image space) data
      
      epi_hyb_roft = mrir_iDFT_freqencode(epi_raw_apod);
      
      
      %==----------------------------------------------------------------==%
      %%% {4}: phase correction

      % assumes each repetition collects new phascor navigators (siemens EPI only)
      linear_fit_coeff = mrir_artifact_ghost_compute(phascor1d_grid);
      epi_hyb_corr = mrir_artifact_ghost_correct(epi_hyb_roft, linear_fit_coeff);
      
      % collapse multiple segments
      if ( FLAG__FIRST_CALL ),
	% probably only need rigorous checks on the first image
	epi_hyb_coll = mrir_multishot_segment_collapse(epi_hyb_corr, meas.evp);
      else,
	epi_hyb_coll = sum(epi_hyb_corr, 8);
      end;
	
      
      %==----------------------------------------------------------------==%
      %%% {5}: 1-D apodization filter to reduce ringing due to truncation

      epi_hyb_apod = mrir_filter_raw_apodize_1d(epi_hyb_coll, 2, 0.15);

      %==----------------------------------------------------------------==%
      %%% {6}: partial fourier processing along phase encoding direction

      if ( meas.prot.ucPhasePartialFourier < 1 ),
%        epi_img_peft = mrir_partial_fourier(epi_hyb_apod, meas.prot, 'homodyne');
        epi_img_peft = mrir_partial_fourier(epi_hyb_apod, meas.prot, 'zero-fill');
      else,
        epi_img_peft = mrir_iDFT_phasencode(epi_hyb_apod);
      end;


      %==----------------------------------------------------------------==%
      %%% back out raw data and crop image data for sending

      %       1 2 3 4 5 6   7 8 9   0
      epi_ksp(:,:,:,1,1,1,rep,1,1,slc) = mrir_fDFT_freqencode(mrir_fDFT_phasencode(epi_img_peft));

      epi_img(:,:,:,1,1,1,rep,1,1,slc) = mrir_image_crop(epi_img_peft);


      FLAG__FIRST_CALL = 0;
    
    end;
  end;
  
  if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_ksp = mrir_image_slice_deinterleave(epi_ksp);
    epi_img = mrir_image_slice_deinterleave(epi_img);
  end;


  epi = mrir_conventional_2d(epi_ksp);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


