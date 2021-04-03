function [epi] = mrir_epi_GS_RawData(meas)

  %==--------------------------------------------------------------------==%

  dims = size(meas.data(:,:,:,1,1,1,:,:,1,:));
  %dims(2) = dims(2) / meas.prot.ucPhasePartialFourier;
  dims(8) = dims(8) / 2;
  
  epi  = complex(zeros(dims), zeros(dims));
  %==--------------------------------------------------------------------==%
  

  global FLAG__FIRST_CALL;
  FLAG__FIRST_CALL = 1;
  
  for rep = 1:meas.evp.NRepMeas,
      for slc = 1:meas.evp.NSlcMeas,
          
          epi_raw_roft = mrir_iDFT_freqencode(meas.data(:,:,:,1,1,1,rep,:,1,slc));
          
          %%% regridding to correct for ramp sampling
          % skip the regridding step if neither trapezoid or sinusoid
          if ( meas.prot.alRegridMode == 1 ),
              epi_raw_grid = epi_raw_roft;
          else
              if ( FLAG__FIRST_CALL ),
                  prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(epi_raw_roft, 1));
              end
              epi_raw_grid = mrir_regrid_trapezoid(epi_raw_roft, meas.prot, prot_trapezoid); 
          end

          epi(:,:,:,1,1,1,rep,1,1,slc) = double(mrir_fDFT_freqencode(epi_raw_grid));
          
          FLAG__FIRST_CALL = 0;
          
      end;
  end;
  
  if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
      epi = mrir_image_slice_deinterleave(epi);     
  end;

  return;



