function evi_img_ind = mrir_evi_GRAPPA(dat, acs, prot, evp, Nx, Ny, Nz, Nk1, Nk2)
%MRIR_EVI_GRAPPA
%
% evi_ind = mrir_evi_GRAPPA(dat, acs, prot, evp, Nx, Ny, Nz, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/07
% $Id: mrir_evi_GRAPPA.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  G = mrir_array_GRAPPA_2d_kernel(acs, evp.NAFLin, evp.NAFPar, Nx, Ny, Nz, Nk1, Nk2);

  dat_raw_full = [];
  for rep = 1:mrir_ice_dimensions(dat, 'rep'),
    dat_raw_redu_rep = dat(:,:,:,:,:,:,rep,:,:,:,:,:,:,:,:,:);
    
    dat_raw_full_rep = mrir_array_GRAPPA_2d_recon(dat_raw_redu_rep, G, evp.NAFLin, evp.NAFPar, evp.NLinMeas, evp.NParMeas, evp.NFirstLin, evp.NFirstPar);
    dat_raw_full = cat(7, dat_raw_full, dat_raw_full_rep);
  
  end;
    
  if ( prot.ucPhasePartialFourier < 1.0 ),
      
    % assume that only partial fourier acquisitions occur in primary phase
    % encoding direction
    evi_roft = mrir_iDFT_freqencode(dat_raw_full);
    evi_paft = mrir_iDFT_partencode(evi_roft);
    
    evi_img = mrir_partial_fourier(evi_paft, prot, 'homodyne', 'pre');
    
    evi_img_ind = mrir_image_crop(evi_img, prot.flReadoutOSFactor);
    
  else,
    
    evi_img_ind = mrir_conventional_3d(dat_raw_full);
    
  end;
  
%  evi_img_rss = mrir_array_combine_rss(evi_img_ind);

  
  return;
    

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_GRAPPA.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
opt.ReturnStruct = 1;
opt.ReadMultipleRepetitions = 0;
opt.PhascorCollapseSegments = 0;


meas_MID331.evp.NFirstRefLin =  7;
meas_MID331.evp.NFirstRefPar =  9;

meas_MID331 = read_meas_dat('meas_MID331_evi32_4x2_pf68_72_FID13884.dat', opt);
