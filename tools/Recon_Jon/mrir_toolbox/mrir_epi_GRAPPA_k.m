function [epi_k_ind, G] = mrir_epi_GRAPPA_k(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2)
%MRIR_EPI_GRAPPA
%
% epi_ind = mrir_epi_GRAPPA(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/07
% $Id: mrir_epi_GRAPPA.m,v 1.1 2008/11/06 01:42:13 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  epi_k_ind = [];
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
      %% Kawin Hack this!!! Only work for R = 2, but can be easily modify
      if(0)
          dat_raw_full_rep_slc = mrir_array_GRAPPA_2d_recon(dat_raw_redu_rep_slc, G{slc}, evp.NLinMeas, evp.NParMeas, evp.NFirstLin, evp.NFirstPar);
      else
          disp('Kawin: Hack to make it fast!!!!')
          for count = 1:length(G{slc}.kernel)
              w(:,:,count) = permute(G{slc}.kernel{count},[2 1]);
          end
          if slc == 1
              s = size(dat_raw_redu_rep_slc);
              s(2) = s(2)*2;
              dat_raw_full_rep_slc = zeros(s);
          end
          K = MultisliceGRAPPA_old(dat_raw_redu_rep_slc,w,[G{slc}.Nsrcx G{slc}.Nsrcy]);
          dat_raw_full_rep_slc(:,1:2:end,:,:,:,:,:,:,:,:) = dat_raw_redu_rep_slc;
          dat_raw_full_rep_slc(:,2:2:end,:,:,:,:,:,:,:,:) = K;
      end
      %%%%%%%%%%%%%%
      
      
      dat_raw_full_slc = cat(7, dat_raw_full_slc, dat_raw_full_rep_slc);
      
    end;

%     
%     if ( prot.ucPhasePartialFourier < 1.0 ),
%       
%       % assume that only partial fourier acquisitions occur in primary phase
%       % encoding direction
%       epi_roft = mrir_iDFT_freqencode(dat_raw_full_slc);
%       
%       epi_img = mrir_partial_fourier(epi_roft, prot, 'zero-fill', 'pre');
%       
%       epi_img_slc = mrir_image_crop(epi_img, prot.flReadoutOSFactor);
%       
%     else,
%       
%       
%       epi_img_slc = mrir_conventional_2d(dat_raw_full_slc, prot);
%       
%     end;
    %     
      
      
    epi_k_slc = dat_raw_full_slc;

    %epi_img_ind = cat(mrir_DIM_SLC, epi_img_ind, epi_img_slc);
    epi_k_ind = cat(mrir_DIM_SLC, epi_k_ind, epi_k_slc);

    
  end;


  if ( strcmp(prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
      %epi_img_ind = mrir_image_slice_deinterleave(epi_img_ind);
      epi_k_ind = mrir_image_slice_deinterleave(epi_k_ind);
  end;

  
  return;
  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi_GRAPPA.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
