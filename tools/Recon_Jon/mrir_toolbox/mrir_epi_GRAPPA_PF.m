function [epi_img_ind, G] = mrir_epi_GRAPPA_PF(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2,PE_type)
%MRIR_EPI_GRAPPA
%
% epi_ind = mrir_epi_GRAPPA(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/07
% $Id: mrir_epi_GRAPPA.m,v 1.1 2008/11/06 01:42:13 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  epi_img_ind = [];
  for slc = 1:mrir_ice_dimensions(dat, 'slc'),

    fprintf(1, '\n');
    disp(sprintf('==> [%s]: slice %02d', mfilename, slc));
    
    if ( ~iscell(acs) ),
      %             1 2 3 4 5 6 7 8 9   0 1 2 3 4 5 6
      acs_slc = acs(:,:,:,:,:,:,:,:,:,slc,:,:,:,:,:,:);
      
      G{slc} = mrir_array_GRAPPA_2d_kernel(acs_slc, evp.NAFLin, evp.NAFPar, Nx, Ny, 1, Nk1, Nk2);
    else,
      G = acs;
      disp('      using pre-computed kernel');
    end;
    
      
    dat_raw_full_slc = [];
    for rep = 1:mrir_ice_dimensions(dat, 'rep'),
      
      %                          1 2 3 4 5 6   7 8 9   0 1 2 3 4 5 6
      dat_raw_redu_rep_slc = dat(:,:,:,:,:,:,rep,:,:,slc,:,:,:,:,:,:);
      
      dat_raw_full_rep_slc = mrir_array_GRAPPA_2d_recon(dat_raw_redu_rep_slc, G{slc}, evp.NLinMeas, evp.NParMeas, evp.NFirstLin, evp.NFirstPar);
      
      dat_raw_full_slc = cat(7, dat_raw_full_slc, dat_raw_full_rep_slc);
      
    end;
    
    Ncol_full = prot.iNoOfFourierColumns;
    Ncol_part = mrir_ice_dimensions(dat_raw_full_slc, 'col');
    
    if  ( prot.ucPhasePartialFourier < 1.0 )
        H = mrir_partial_fourier(dat_raw_full_slc, prot, 'zero-fill', 'pre');
    else
        H = mrir_iDFT_phasencode(dat_raw_full_slc);
    end

    if  ( Ncol_full ~= Ncol_part )
        lPhaseEncodingLinesBackup = prot.lPhaseEncodingLines;
        prot.lPhaseEncodingLines = prot.iNoOfFourierColumns; %if fully acquired
        ucPhasePartialFourierBackup = prot.ucPhasePartialFourier;
        prot.ucPhasePartialFourier = size(H,1)/prot.iNoOfFourierColumns ;
        
        if PE_type == 1
            %epi_img = permute(mrir_partial_fourier(permute(H,[2 1 3:16]), prot, 'hymodyne', 'pre'),[2 1 3:16]);
            epi_img = permute(mrir_partial_fourier(permute(H,[2 1 3:16]), prot, 'zero-fill', 'pre'),[2 1 3:16]);
        else
            %epi_img = permute(mrir_partial_fourier(permute(H,[2 1 3:16]), prot, 'hymodyne', 'post'),[2 1 3:16]);
            epi_img = permute(mrir_partial_fourier(permute(H,[2 1 3:16]), prot, 'zero-fill', 'pre'),[2 1 3:16]);
        end
        
        prot.lPhaseEncodingLines = lPhaseEncodingLinesBackup;
        prot.ucPhasePartialFourier = ucPhasePartialFourierBackup; 
    else
        epi_img = mrir_iDFT_freqencode(H);
    end

    epi_img_slc = mrir_image_crop(epi_img, prot.flReadoutOSFactor);

    epi_img_ind = cat(mrir_DIM_SLC, epi_img_ind, epi_img_slc);

    
  end;


  if ( strcmp(prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_img_ind = mrir_image_slice_deinterleave(epi_img_ind);
  end;

  
  return;
  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi_GRAPPA.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
