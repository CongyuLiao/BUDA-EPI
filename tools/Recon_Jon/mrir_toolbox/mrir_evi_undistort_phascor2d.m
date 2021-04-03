function [acs_raw_undistort, acs_img_undistort] = mrir_evi_undistort_phascor2d(meas)
%MRIR_EVI_UNDISTORT_PHASCOR2D
%
% acs_undistort = mrir_evi_undistort_phascor2d(meas)

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
      [ repmat((1:R1)', R2,1), reshape(repmat((1:R2), R1,1), [], 1) ];


  acs_raw = meas.patrefscan;
  nav_raw = meas.data_phascor2d;

  % apply conventional 1D phase correction to reduce ghosting before
  % computing 2D correction

  phascor1d_coeff = mrir_artifact_ghost_compute(meas.patrefscan_phascor);
  prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(acs_raw, 1));

  %===----

  % transform into hybrid space to apply 1D correction
  acs_hyb_roft = mrir_iDFT_freqencode(meas.patrefscan);
  acs_hyb_corr = mrir_artifact_ghost_correct(acs_hyb_roft, phascor1d_coeff);

  acs_hyb_coll = mrir_multishot_segment_collapse(acs_hyb_corr, meas.evp);
%  acs_hyb_grid = mrir_regrid_trapezoid(acs_hyb_coll, meas.prot, prot_trapezoid, 1);
 
  %===----

  acs_hyb_peft = mrir_iDFT_phasencode(acs_hyb_coll);   % *** not used ***
  acs_img_paft = mrir_iDFT_partencode(acs_hyb_peft);   % *** not used ***
  acs_rss = mrir_array_combine_rss(acs_img_paft);      % *** not used ***

  
  %===----

  % transform into hybrid space to apply 1D correction
  nav_hyb_roft = mrir_iDFT_freqencode(meas.data_phascor2d);
  nav_hyb_corr = mrir_artifact_ghost_correct(nav_hyb_roft, phascor1d_coeff);

  nav_hyb_coll = sum(nav_hyb_corr, mrir_DIM_SEG);
%  nav_hyb_grid = mrir_regrid_trapezoid(nav_hyb_coll, meas.prot, prot_trapezoid, 1);


  %===----

  for ind = 1:Nshots,

    visits = ( (ind-1)*Nvisit_per_shot+1 ) : ( (ind-0)*Nvisit_per_shot+0 );

    % calculate k-space lines and partitions in fully-sampled data that
    % correspond in time to the navigators collected in this shot (based the
    % pre-calculated first line and first partition for each shot).
    lins = patrefscan_FirstLin_FirstPar_shots(ind, 1):R1:Nlin;
    pars = patrefscan_FirstLin_FirstPar_shots(ind, 2):R2:Npar;

    acs_hyb_extract = acs_hyb_coll(:,lins,:,:,:,:,:,:,pars);
    nav_hyb_extract = nav_hyb_coll(:,lins,:,:,:,:,:,:,visits);

    acs_hyb_extract_peft = mrir_iDFT_phasencode(acs_hyb_extract);
    nav_hyb_extract_peft = mrir_iDFT_phasencode(nav_hyb_extract);
    
    % pick first partition's worth of navigators as reference
    phase_ref = angle(nav_hyb_extract_peft(:,:,:,:,:,:,:,:,1,:,:,:,:,:,:,:));

    % since "true" phase from image projection is assumed constant over
    % time, any phase differences are from off-resonance dephasing, so
    % calculate the phase difference between partitions and use this as
    % phase adjustment.
    phase_adj = angle(nav_hyb_extract_peft) ...
        - repmat(phase_ref, [1,1,1,1,1,1,1,1,Nvisit_per_shot,1,1,1,1,1,1,1]);

    % this correction, when applied, should give each partition the phase of
    % the reference partition
    cplx_corr = exp(-i*phase_adj);
    
%    disp('flipping...'); cplx_corr = flipdim(cplx_corr, 2);

    acs_hyb_extract_undistort = acs_hyb_extract_peft .* cplx_corr;
    nav_hyb_extract_undistort = nav_hyb_extract_peft .* cplx_corr;

    acs_hyb_extract_undistort_roft = mrir_fDFT_phasencode(acs_hyb_extract_undistort);
    nav_hyb_extract_undistort_roft = mrir_fDFT_phasencode(nav_hyb_extract_undistort);
    
    
    % stuff corrected data back into appropriate lines and partitions
    acs_hyb_undistort_roft(:,lins,:,:,:,:,:,:,pars)   = acs_hyb_extract_undistort_roft;
    nav_hyb_undistort_roft(:,lins,:,:,:,:,:,:,visits) = nav_hyb_extract_undistort_roft;

    %==------------------------------------------------------------------==%
    if ( DEBUG ),

      cha = 02;
      
      % generate some debugging figures
      mrir_evi_undistort_phascor2d__DEBUG(nav_hyb_extract_peft, nav_hyb_extract_undistort, acs_hyb_extract_peft, acs_hyb_extract_undistort, lins, Nvisit_per_shot, cha, ind);
      
    end;
    %==------------------------------------------------------------------==%

  end;  %% FOR loop



  %===----

  acs_raw_undistort = mrir_fDFT_freqencode(acs_hyb_undistort_roft);

  acs_hyb_undistort_grid = mrir_regrid_trapezoid(acs_hyb_undistort_roft, meas.prot, prot_trapezoid, 1);

  acs_hyb_undistort_peft = mrir_iDFT_phasencode(acs_hyb_undistort_grid);
  acs_img_undistort_paft = mrir_iDFT_partencode(acs_hyb_undistort_peft);
  
  acs_img_undistort = mrir_image_crop(acs_img_undistort_paft);
  
  acs_rss_undistort = mrir_array_combine_rss(acs_img_undistort);


  %===----

  nav_hyb_undistort = mrir_fDFT_freqencode(nav_hyb_undistort_roft);

  nav_hyb_undistort_grid = mrir_regrid_trapezoid(nav_hyb_undistort_roft, meas.prot, prot_trapezoid, 1);

  nav_hyb_undistort_peft = mrir_iDFT_phasencode(nav_hyb_undistort_grid);
  nav_img_undistort_paft = mrir_iDFT_partencode(nav_hyb_undistort_peft);
  
  nav_img_undistort = mrir_image_crop(nav_img_undistort_paft);
  
  nav_rss_undistort = mrir_array_combine_rss(nav_img_undistort);

  
  
  
  return;



%**************************************************************************%
function h = mrir_evi_undistort_phascor2d__DEBUG(nav_hyb_extract_peft, nav_hyb_extract_undistort, acs_hyb_extract_peft, acs_hyb_extract_undistort, lins, N, cha, ind)

  h(1) = figure('name', sprintf('%s (shot %03d)', mfilename, ind));

  %      subplot(3,2,2);
  %      imagesc(mod(mrir_image_mosaic(squeeze(phase_adj(:,:,cha,1,1,1,1,1,:)), [1, N]), 2*pi))
  %      set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on');
  %      axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  
  subplot(2,2,1);
  imagesc(angle(mrir_image_mosaic(squeeze(nav_hyb_extract_peft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  %      title('"phascor2d" navigators -- original')
  
  
  subplot(2,2,2);
  imagesc(angle(mrir_image_mosaic(squeeze(nav_hyb_extract_undistort(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  %      title('"phascor2d" navigators -- corrected')
  
  
  subplot(2,2,3);
  imagesc(angle(mrir_image_mosaic(squeeze(acs_hyb_extract_peft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  title('"patrefscan" ACS lines -- original')
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  
  
  subplot(2,2,4);
  imagesc(angle(mrir_image_mosaic(squeeze(acs_hyb_extract_undistort(:,:,cha,1,1,1,1,1,:)), [1, N])));
  title('"patrefscan" ACS lines -- corrected')
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  
  set(suptitle(sprintf('%s (shot %03d)', mfilename, ind)), 'Interpreter', 'none');

  
  %------------------------------------------------------------------------%

  h(2) = figure('name', sprintf('%s (shot %03d)', mfilename, ind));
  
  nav_hyb_extract_paft = mrir_iDFT_partencode(nav_hyb_extract_peft);
  nav_hyb_extract_undistort_paft = mrir_iDFT_partencode(nav_hyb_extract_undistort);

  acs_hyb_extract_paft = mrir_iDFT_partencode(acs_hyb_extract_peft);
  acs_hyb_extract_undistort_paft = mrir_iDFT_partencode(acs_hyb_extract_undistort);
  
  subplot(2,2,1);
  imagesc(abs(mrir_image_mosaic(squeeze(nav_hyb_extract_paft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  %      title('"phascor2d" navigators -- original')

  subplot(2,2,2);
  imagesc(abs(mrir_image_mosaic(squeeze(nav_hyb_extract_undistort_paft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  %      title('"phascor2d" navigators -- corrected')
 
  subplot(2,2,3);
  imagesc(abs(mrir_image_mosaic(squeeze(acs_hyb_extract_paft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  title('"phascor2d" navigators -- original')

  subplot(2,2,4);
  imagesc(abs(mrir_image_mosaic(squeeze(acs_hyb_extract_undistort_paft(:,:,cha,1,1,1,1,1,:)), [1, N])));
  set(gca, 'XTick', length(lins):length(lins):(length(lins)*N), 'XGrid', 'on', 'GridLineStyle', '-');
  axis image; xlabel('phase encoding'); ylabel('frequency encoding');
  title('"phascor2d" navigators -- corrected')
  
  set(suptitle(sprintf('%s (shot %03d)', mfilename, ind)), 'Interpreter', 'none');
  
  
  
  
  
  
  
  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_distort_B0_phaseref.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
