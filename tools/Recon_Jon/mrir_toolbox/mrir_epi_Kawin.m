function [epi, epi_raw] = mrir_epi_K(meas)
%MRIR_EPI
%
% epi = mrir_epi(meas_struct)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/16
% $Id: mrir_epi.m,v 1.5 2007/11/17 00:47:38 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims = size(meas.data(:,:,:,1,1,1,:,:,1,:));
  dims(1) = dims(1) / 2;
  dims(2) = dims(2) / meas.prot.ucPhasePartialFourier;
  dims(8) = dims(8) / 2;
  
  epi     = complex(zeros(dims), zeros(dims));
  

  %==--------------------------------------------------------------------==%

  global FLAG__FIRST_CALL;
  FLAG__FIRST_CALL = 1;

  for rep = 1:meas.evp.NRepMeas,
    for slc = 1:meas.evp.NSlcMeas,

      if ( meas.evp.NRepMeas > 1 ),
	disp(sprintf('  rep %03d', rep));
      end;
      
      epi_raw_roft = mrir_iDFT_freqencode(meas.data(:,:,:,1,1,1,rep,:,1,slc));

      %==----------------------------------------------------------------==%
      %%% phase correction
      
      % assumes each repetition collects new phascor navigators (siemens EPI only)
      linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d(:,:,:,1,1,1,rep,:,1,slc));
      epi_raw_corr = mrir_artifact_ghost_correct(epi_raw_roft, linear_fit_coeff);

% for kawin
% epi_raw_corr = epi_raw_roft;

      if ( FLAG__FIRST_CALL ),
	% probably only need rigorous checks on the first image
	epi_raw_corr = mrir_multishot_segment_collapse(epi_raw_corr, meas.evp);
      else,
	epi_raw_corr = sum(epi_raw_corr, 8);
      end;


      %==----------------------------------------------------------------==%
      %%% regridding to correct for ramp sampling

      % skip the regridding step if neither trapezoid or sinusoid
      if ( meas.prot.alRegridMode == 1 ),
	epi_raw_grid = epi_raw_corr;
      else,

	if ( FLAG__FIRST_CALL ),
	  prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(epi_raw_corr, 1));
	end;

	% regrid *after* phase correction
	epi_raw_grid = mrir_regrid_trapezoid(epi_raw_corr, meas.prot, prot_trapezoid);

      end;

      %==----------------------------------------------------------------==%
      %%% partial fourier processing along phase encoding direction

      if ( meas.prot.ucPhasePartialFourier < 1 ),
        epi_img_peft = mrir_partial_fourier(epi_raw_grid, meas.prot, 'zero-fill', 'pre');
      else,
        epi_img_peft = mrir_iDFT_phasencode(epi_raw_grid);
      end;


      %==----------------------------------------------------------------==%
      %%% back out raw data and crop image data for sending

       %       1 2 3 4 5 6   7 8 9   0
      %
      epi_raw(:,:,:,1,1,1,rep,1,1,slc) = mrir_fDFT_freqencode(mrir_fDFT_phasencode(epi_img_peft));

      epi_img(:,:,:,1,1,1,rep,1,1,slc) = mrir_image_crop(epi_img_peft);


% Kawin Setsompop, May 14 2009, return double percision result
      %epi_raw(:,:,:,1,1,1,rep,1,1,slc) = double(mrir_fDFT_freqencode(mrir_fDFT_phasencode(epi_img_peft)));
      %epi_img(:,:,:,1,1,1,rep,1,1,slc) = double(mrir_image_crop(epi_img_peft));

      

      FLAG__FIRST_CALL = 0;

    end;
  end;

  if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    %    epi_raw = mrir_image_slice_deinterleave(epi_raw);
    %    epi_img = mrir_image_slice_deinterleave(epi_img);
    epi_img = mrir_image_slice_deinterleave(epi_img);
  end;


%  epi = mrir_conventional_2d(epi_raw);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


