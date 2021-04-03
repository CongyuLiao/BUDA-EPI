function [epi_img_ind, G] = mrir_epi_GRAPPA_v3(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2, apodization_para,eta,chi,col)
%MRIR_EPI_GRAPPA
%
% epi_ind = mrir_epi_GRAPPA(dat, acs, prot, evp, Nx, Ny, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/07
% $Id: mrir_epi_GRAPPA.m,v 1.1 2008/11/06 01:42:13 jonnyreb Exp $


% Kawin Setsompop, April 24 2010:  modify to add apodization at the end
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
      
      %G{slc} = mrir_array_GRAPPA_2d_kernel(acs_slc, evp.NAFLin, evp.NAFPar, Nx, Ny, 1, Nk1, Nk2);
      %G{slc} = mrir_array_GRAPPA_2d_kernel_improved_dev(acs_slc, evp.NAFLin, evp.NAFPar, Nx, Ny, 1, Nk1, Nk2);
      G{slc} = mrir_array_GRAPPA_2d_kernel_improved_dev(acs_slc, evp.NAFLin, evp.NAFPar, Nx, Ny, 1, Nk1, Nk2,eta,chi,col);     
    else,
      G = acs;
      disp('      using pre-computed kernel');
    end;
    
    
    dat_raw_full_slc = [];
    NCol = mrir_ice_dimensions(dat, 'col');
    NCha = mrir_ice_dimensions(dat, 'cha');
    for rep = 1:mrir_ice_dimensions(dat, 'rep'),
      
      %                          1 2 3 4 5 6   7 8 9   0 1 2 3 4 5 6
      dat_raw_redu_rep_slc = dat(:,:,:,:,:,:,rep,:,:,slc,:,:,:,:,:,:);
      METHOD = 'conv2';
      switch METHOD
          case 'imagedomain'
              tic
              [i_conv, dat_raw_full_rep_slc] = mrir_array_GRAPPA_recon_imagedomain(dat_raw_redu_rep_slc, G{slc}, evp);     
              toc
          case 'old'
              dat_raw_full_rep_slc = mrir_array_GRAPPA_2d_recon(dat_raw_redu_rep_slc, G{slc}, evp.NLinMeas, evp.NParMeas, evp.NFirstLin, evp.NFirstPar);
          case 'conv2'

              %[i_conv, k_conv] = mrir_array_GRAPPA_recon_kspace_convolve(dat_raw_redu_rep_slc, G{slc}, evp);
              %[kernel_conv, kernel_corr] = mrir_array_GRAPPA_conv_kernel(G{slc});
              
              tic
              % its AMAZING that this works...
              k_datlines1 = evp.NFirstLin : G{slc}.R1 : evp.NLinMeas;    
              k_datlines2 = evp.NFirstPar : G{slc}.R2 : evp.NParMeas;
              
              k_fill = complex(zeros( NCol, evp.NLinMeas, NCha, 1, 1, 1, 1, 1, evp.NParMeas, class(dat_raw_redu_rep_slc) ));
              k_fill(:, k_datlines1, :, 1,1,1,1,1, k_datlines2) = dat_raw_redu_rep_slc;
              
              k = mrir_array_GRAPPA_conv_kernel_v2(G{slc}); % not quite correct as shift the resulting k-space by one....
              
              k_conv_cha = complex(zeros(NCol, evp.NLinMeas, NCha, NCha));
              
              for trg = 1:NCha,
                  for src = 1:NCha,
                      k_conv_cha(:,:,src,trg) = conv2( k_fill(:,:,src), k(:,:,1,src,trg), 'same');
                  end;
              end;
              
              dat_raw_full_rep_slc = squeeze(sum(k_conv_cha, 3));
              toc

          case 'conv2_fast'
              % still need to check generalization to more than R_inplane =
              % 2 -> check if the "w = G{slc}.kernel{1}(ind:Ry:end,:).'" is
              % correct and cases where Nlin is not divisible by R_inplane
              tic
              
              if slc == 1 && rep == 1
                  dat_raw_full_rep_slc = complex(zeros( NCol, evp.NLinMeas, NCha, 1, 1, 1, 1, 1, evp.NParMeas, class(dat_raw_redu_rep_slc) ));
                  R1 = G{1}.R1; Ry = max([1, (R1-1)]);
              end
              dat_raw_full_rep_slc(:,evp.NFirstLin:R1:end,:,:,:,:,:,:,:,:) = dat_raw_redu_rep_slc;
              %dat_raw_full_rep_slc(:,1:R1:end,:,:,:,:,:,:,:,:) = dat_raw_redu_rep_slc;
              
              ind = 0;
              K_current = complex(zeros(size(dat_raw_redu_rep_slc)));
              for ind_y = 1:Ry,                  
                  ind = ind + 1;
                  w = G{slc}.kernel{1}(ind:Ry:end,:).';   % G:  [NCha * (R1-1) * (R2-1)] x [NCha * Nsrcx * Nsrcy * Nsrcz
                  for ChCount = 1:NCha
                      Wcurrent = reshape(w(:,ChCount),G{slc}.Nsrcx,G{slc}.Nsrcy,NCha);
                      for ChCount2 = 1:NCha
                          WcurrentCh2 = Wcurrent(end:-1:1,end:-1:1,ChCount2);
                          K_current(:,:,ChCount,:,:,:,:,:,:,:) =  K_current(:,:,ChCount,:,:,:,:,:,:,:) + conv2(dat_raw_redu_rep_slc(:,:,ChCount2,:,:,:,:,:,:,:),WcurrentCh2,'same') ;
                      end
                  end
                  
                  startIndex = ind+evp.NFirstLin;
                  if startIndex > R1
                      startIndex = rem(startIndex,R1);
                      %startIndex2 = 
                  end
                  % startIndex = ind+1;
                  NlinCurrent = size(dat_raw_full_rep_slc(:,startIndex:R1:end,:,:,:,:,:,:,:,:),2);
                  %dat_raw_full_rep_slc(:,startIndex:R1:end,:,:,:,:,:,:,:,:) = K_current(:,1:NlinCurrent,:,:,:,:,:,:,:,:);
                  dat_raw_full_rep_slc(:,startIndex:R1:end,:,:,:,:,:,:,:,:) = K_current(:,1:NlinCurrent,:,:,:,:,:,:,:,:);
              end;  
              toc
      end
      %%%%%%%%%%%%%%
      
      
      dat_raw_full_slc = cat(7, dat_raw_full_slc, dat_raw_full_rep_slc);
      
    end;

    
    if ( prot.ucPhasePartialFourier < 1.0 ),
        
        % assume that only partial fourier acquisitions occur in primary phase
        % encoding direction
        epi_roft = mrir_iDFT_freqencode(dat_raw_full_slc);
        
        epi_img = mrir_partial_fourier(epi_roft, prot, 'zero-fill', 'pre');
        
        
        if (apodization_para ~= 0) % apotization!!!
            epi_k = mrir_fDFT_freqencode(mrir_fDFT_phasencode(epi_img));
            epi_k_apod = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(epi_k, 1, apodization_para),  2, apodization_para) ;
            epi_img_slc = mrir_conventional_2d(epi_k_apod,prot);
        else
            epi_img_slc = mrir_image_crop(epi_img, prot.flReadoutOSFactor);
        end
    else,
        if (apodization_para ~= 0) % apotization!!!
            dat_raw_full_slc = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(dat_raw_full_slc, 1, apodization_para),  2, apodization_para) ;
        end
        epi_img_slc = mrir_conventional_2d(dat_raw_full_slc, prot);
    end;
    
    
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
