function evi_undistort = mrir_evi_undistort_phascor2d(meas)
%MRIR_EVI_UNDISTORT_PHASCOR2D
%

% note that this function contains logic that is specific to the MGH EVI
% acquisition specification and will not generalize smoothly to other
% schemes. thus, any updates to the EVI specification must be incorporated
% here as well.

% van der Zwaag W, Francis S, Bowtell R. "Improved echo volumar imaging
% (EVI) for functional MRI". Magn Reson Med. 2006 Dec;56(6):1320-7. PMID:
% 17089364.


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2008/may/14
% $Id: mrir_evi_distort_B0_phaseref.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  R1 = meas.evp.NAFLin;
  R2 = meas.evp.NAFPar;

  % "acquisitions" is same as number of shots
  Nshots = R1 * R2;

  % numbers of lines and paritions in the ACS scans
  Nlin = meas.evp.NRefLin;
  Npar = meas.evp.NRefPar;

  % each partition is visited R1 times during phascor2d navigators
  Nvisit_per_par = R1;

  Nvisit_per_shot = Npar * Nvisit_per_par / Nshots;


  % since partition counter in the phascor2d navigators is independent of
  % the partitions in k-space of the corresponding ACS or image scans, to
  % interpret which navigators correspond in time to which lines and
  % partitions here we pre-calculate only the first line and first partition
  % visited in each shot of the ACS acquisition and later calculate all of
  % the lines and partitions to interpret the navigators. this relies
  % heavily on the known acquisition scheme of the MGH EVI pulse sequence
  % (e.g., that the ACS lines always start with line 0 and partition 0), but
  % it may be possible to generalize.

  patrefscan_FirstLin_FirstPar_shots = ...
      [repmat((1:R1)', R2,1), reshape(repmat((1:R2), R1,1), [], 1) ]

  phascor1d_coeff = mrir_artifact_ghost_compute(meas.patrefscan_phascor);


  %===----

  acs_hyb_roft = mrir_iDFT_freqencode(meas.patrefscan);
  acs_hyb_corr = mrir_artifact_ghost_correct(acs_hyb_roft, phascor1d_coeff);
  acs_hyb_coll = mrir_multishot_segment_collapse(acs_hyb_corr, meas.evp);

  acs_hyb_peft = mrir_iDFT_phasencode(acs_hyb_coll);


  %===----

  acs_img_paft = mrir_iDFT_partencode(acs_hyb_peft);
  acs_img = mrir_image_crop(acs_img_paft);
  acs_rss = mrir_array_combine_rss(acs_img);


  %===----

  nav_hyb_roft = mrir_iDFT_freqencode(meas.data_phascor2d);
  nav_hyb_corr = mrir_artifact_ghost_correct(nav_hyb_roft, phascor1d_coeff);



  for ind = 1:Nshots,

    visits = ( (ind-1)*Nvisit_per_shot+1 ) : ( (ind-0)*Nvisit_per_shot+0 );

    % calculate k-space lines and partitions in fully-sampled data that
    % correspond in time to the navigators collected in this shot (based the
    % pre-calculated first line and first partition for each shot).
    lins = patrefscan_FirstLin_FirstPar_shots(ind, 1):R1:Nlin;
    pars = patrefscan_FirstLin_FirstPar_shots(ind, 2):R2:Npar;

    acs_hyb_peft_extract = acs_hyb_peft(:,lins,:,:,:,:,:,:,pars);

    nav_hyb_corr_extract = nav_hyb_corr(:,lins,:,:,:,:,:,:,visits);
    nav_hyb_coll_extract = sum(nav_hyb_corr_extract, mrir_DIM_SEG);
    nav_hyb_peft_extract = mrir_iDFT_phasencode(nav_hyb_coll_extract);

    phase_ref = angle(nav_hyb_peft_extract(:,:,:,:,:,:,:,:,1,:,:,:,:,:,:,:));
    phase_adj = angle(nav_hyb_peft_extract) ...
        - repmat(phase_ref, [1,1,1,1,1,1,1,1,Nvisit_per_shot,1,1,1,1,1,1,1]);

    cplx_corr = exp(-i*phase_adj);


    nav = nav_hyb_peft_extract .* cplx_corr;
    acs = acs_hyb_peft_extract .* cplx_corr;

    nav_hyb_peft_undistort(:,lins,:,:,:,:,:,:,visits) = nav;
    acs_hyb_peft_undistort(:,lins,:,:,:,:,:,:,pars)   = acs;


    %===----

    if ( DEBUG ),

      figure('name', sprintf('%s  shot %03d', mfilename, ind));
      subplot(2,2,1);
      imagesc(angle(mrir_image_mosaic(squeeze(nav_hyb_peft_extract(:,:,1,1,1,1,1,1,:)), [1, 6])));
      axis image;
      title('"phascor2d" navigators')

      subplot(2,2,2);
      imagesc(angle(mrir_image_mosaic(squeeze(nav(:,:,1,1,1,1,1,1,:)), [1, 6])));
      axis image;
      title('"phascor2d" navigators -- corrected')

      subplot(2,2,3);
      imagesc(angle(mrir_image_mosaic(squeeze(acs_hyb_peft_extract(:,:,1,1,1,1,1,1,:)), [1, 6])));
      title('"patrefscan" ACS lines')
      axis image;

      subplot(2,2,4);
      imagesc(angle(mrir_image_mosaic(squeeze(acs(:,:,1,1,1,1,1,1,:)), [1, 6])));
      title('"patrefscan" ACS lines -- corrected')
      axis image;

      set(suptitle(sprintf('%s  shot %03d', mfilename, ind)), 'Interpreter', 'none');

    end;


  end;


  %  nav_img_paft = mrir_iDFT_partencode(nav_hyb_peft);
  %  nav_img = mrir_image_crop(nav_img_paft);
  %  nav_rss = mrir_array_combine_rss(nav_img);

  acs_img_paft_undistort = mrir_iDFT_partencode(acs_hyb_peft_undistort);
  acs_img_undistort = mrir_image_crop(acs_img_paft_undistort);
  acs_rss_undistort = mrir_array_combine_rss(acs_img_undistort);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_distort_B0_phaseref.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


  %    acs       = meas.patrefscan(    :,lins,:,:,:,:,:,:,pars);

    %                               1    2 3 4 5 6 7 8       9 0 1 2 3 4 5 6
%    phascor2d = meas.data_phascor2d(:,lins,:,:,:,:,:,:,visits,:,:,:,:,:,:,:);

    %                           1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
%    phase_ref = angle(phascor2d(:,:,:,:,:,:,:,:,1,:,:,:,:,:,:,:));


    %                                                  1 2 3 4 5 6 7 8            9 0 1 2 3 4 5 6
%    phase_corr = angle(phascor2d) - repmat(phase_ref, [1,1,1,1,1,1,1,1,Nvisit_per_shot,1,1,1,1,1,1,1]);


    %%%    phase_diff = diff(angle(phascor2d), 1, mrir_DIM_PAR);
    %%%    phase_corr = cat(mrir_DIM_PAR, zeros(size(sum(phascor2d,mrir_DIM_PAR))), phase_diff);

