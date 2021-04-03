function [evi, evi_ksp, evi_img] = mrir_evi(meas)
%MRIR_EVI
%
% evi_img = mrir_evi(meas_struct)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/16
% $Id: mrir_evi.m,v 1.3 2007/11/17 00:47:25 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  global FLAG__FIRST_CALL;
  FLAG__FIRST_CALL = 1;

  % initialize variables that accumulate across partitions
  evi_ksp = [];
  evi_img = [];

  for rep = 1:meas.evp.NRepMeas,

    % extract data      1 2 3 4 5 6   7 8 9 0 1 2 3 4 5 6
    evi_raw = meas.data(:,:,:,1,1,1,rep,:,:,1,1,1,1,1,1,1);

    % extract navigators for 1D phase correction  7 8 9 0 1 2 3 4 5 6
    phascor1d = meas.data_phascor1d(:,:,:,1,1,1,rep,:,:,1,1,1,1,1,1,1);


    %==----------------------------------------------------------------==%

    evi_raw_apod = mrir_filter_raw_apodize_1d(evi_raw, mrir_DIM_COL, 0.15);
    evi_hyb_roft = mrir_iDFT_freqencode(evi_raw_apod);


    %==----------------------------------------------------------------==%
    %%% phase correction

    if ( meas.evp.NSegMeas > 2 ),
      disp(sprintf('==> [%s]:  multi-shot data detected!', mfilename));
    end;

    % 2008/apr/12: assuming that if data was collected with N shots, that
    % phase correction navigators also collected with N shots. this
    % information may be redundant, however, so not clear that we need to go
    % to this trouble.
    for ind_shot = 1:(meas.evp.NSegMeas/2),

      % calculate segment index range for forward and reverse lines in
      % this shot
      seg_fwd = ind_shot*2 - 1;
      seg_rev = ind_shot*2 - 0;

      %                          1 2 3 4 5 6 7               8 9 0 1 2 3 4 5 6
      phascor1d_shot = phascor1d(:,:,:,:,:,:,:,seg_fwd:seg_rev,:,1,1,1,1,1,1,1);
      %                                1 2 3 4 5 6 7               8 9 0 1 2 3 4 5 6
      evi_hyb_roft_shot = evi_hyb_roft(:,:,:,:,:,:,:,seg_fwd:seg_rev,:,1,1,1,1,1,1,1);


      % assumes each repetition collects new phascor navigators (siemens EVI only)
      linear_fit_coeff = mrir_artifact_ghost_compute(phascor1d_shot);
      evi_hyb_corr_shot = mrir_artifact_ghost_correct(evi_hyb_roft_shot, linear_fit_coeff);

      evi_hyb_corr(:,:,:,:,:,:,:,seg_fwd:seg_rev,:,1,1,1,1,1,1,1) = evi_hyb_corr_shot;

    end;


    if ( FLAG__FIRST_CALL ),
      % probably only need rigorous checks on the first image
      evi_hyb_coll = mrir_multishot_segment_collapse(evi_hyb_corr, meas.evp);
    else,
      evi_hyb_coll = sum(evi_hyb_corr, 8);
    end;


    %==----------------------------------------------------------------==%
    %%% regridding to correct for ramp sampling

    % skip the regridding step if neither trapezoid or sinusoid
    if ( meas.prot.alRegridMode == 1 ),
      evi_hyb_grid = evi_hyb_coll;
    else,

      if ( FLAG__FIRST_CALL ),
        prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(evi_raw, 1));
      end;

      % regrid *after* phase correction
      evi_hyb_grid = mrir_regrid_trapezoid(evi_hyb_coll, meas.prot, prot_trapezoid, 1);

    end;


    %==----------------------------------------------------------------==%
    %%% partial fourier processing along primary phase encoding direction

    evi_hyb_apod = mrir_filter_raw_apodize_1d(evi_hyb_grid, mrir_DIM_LIN, 0.15);

    if ( meas.prot.ucPhasePartialFourier < 1 ),
      evi_hyb_peft = mrir_partial_fourier(evi_hyb_apod, meas.prot, 'homodyne');
    else,
      evi_hyb_peft = mrir_iDFT_phasencode(evi_hyb_apod);
    end;

    evi_img_apod = mrir_filter_raw_apodize_1d(evi_hyb_peft, mrir_DIM_PAR, 0.0);


    %==----------------------------------------------------------------==%
    %%% back out raw data and crop image data for sending

    evi_ksp = cat(mrir_DIM_REP, evi_ksp, ...
                  mrir_fDFT_freqencode(mrir_fDFT_phasencode(evi_img_apod)));

    % assumes no partial fourier in secondary phase (or partition)
    % encoding direction!
    evi_img_paft = mrir_iDFT_partencode(evi_img_apod);

    evi_img = cat(mrir_DIM_REP, evi_img, ...
                  mrir_image_crop(evi_img_paft));

    FLAG__FIRST_CALL = 0;

  end;


  evi = mrir_conventional_3d(evi_ksp);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


