function [epi, epi_raw, epi_img] = mrir_epi_simple(meas)
%MRIR_EPI_SIMPLE
%
% epi_img = mrir_epi_simple(meas_struct)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/16
% $Id: mrir_epi_simple.m,v 1.1 2007/11/16 23:50:38 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  NLin = meas.evp.NLinMeas;
  NCol = meas.evp.NColMeas;
  
  for rep = 1:meas.evp.NRepMeas,
    for slc = 1:meas.evp.NSlcMeas,

      epi_raw_roft = mrir_iDFT_freqencode(meas.data(:,:,:,1,1,1,rep,:,1,slc));

      %==----------------------------------------------------------------==%
      %%% phase correction
      
      % assumes each repetition collects new phascor navigators (siemens EPI only)
      linear_fit_coeff = mrir_artifact_ghost_compute(sum(meas.data_phascor1d(:,:,:,1,1,1,rep,:,1,slc),8));
      epi_raw_corr = mrir_multishot_segment_collapse(mrir_artifact_ghost_correct(epi_raw_roft, linear_fit_coeff));


      %==----------------------------------------------------------------==%
      %%% partial fourier processing

      if ( meas.prot.ucPhasePartialFourier < 1 ),
	epi_img_peft = mrir_partial_fourier(epi_raw_corr, meas.prot, 'homodyne');
      else,
	epi_img_peft = mrir_iDFT_phasencode(epi_raw_corr);
      end;
      
      %==----------------------------------------------------------------==%
      %%% regridding to correct for ramp sampling

      % regrid *after* phase correction
      epi_img_grid = mrir_regrid_trapezoid(epi_img_peft, meas.prot);

      %       1 2 3 4 5 6   7 8 9   0
      epi_raw(:,:,:,1,1,1,rep,1,1,slc) = mrir_fDFT_freqencode(mrir_fDFT_phasencode(epi_img_grid));

      epi_img(:,:,:,1,1,1,rep,1,1,slc) = mrir_image_crop(epi_img_grid);

    end;
  end;

  if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_raw = mrir_image_slice_deinterleave(epi_raw);
    epi_img = mrir_image_slice_deinterleave(epi_img);
  end;


  epi = mrir_conventional_2d(epi_raw);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi_simple.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


