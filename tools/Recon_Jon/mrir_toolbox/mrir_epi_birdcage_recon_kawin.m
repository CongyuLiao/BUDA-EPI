function [epi_rss, epi_birdcage_01] = mrir_epi_birdcage_recon_kawin(filename, channel_phases)
%MRIR_EPI_BIRDCAGE_RECON_KAWIN  reconstruct segmented EPI and synthesize 1st birdcage mode
%
% [epi_rss, epi_birdcage_01] = mrir_epi_birdcage_recon_kawin(filename, channel_phases)
%
% "channel_phases" are assumed ordered using Kawin's convention, and are
% automatically resorted to match the channels. also, it is assumed that the
% channels appear in order in the raw data file.

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/23
% $Id: mrir_epi_birdcage_recon_kawin.m,v 1.1 2007/10/23 18:44:36 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  % read in data

  opt.ReturnStruct = 1;

  % cannot reorder channels on 7T
  opt.CanonicalReorderCoilChannels = 0;

  % read data and protocol into "meas" struct
  meas = read_meas_dat(filename, opt);


  NCol = mrir_ice_dimensions(meas.data, 'col');
  NLin = mrir_ice_dimensions(meas.data, 'lin');
  NCha = mrir_ice_dimensions(meas.data, 'cha');
  NSeg = mrir_ice_dimensions(meas.data, 'seg');


  %==--------------------------------------------------------------------==%
  % standard EPI recon for segmented acquisition and phase correction for
  % ghosting

  epi_raw_roft = mrir_iDFT_freqencode(meas.data);

  lines_per_segment = NLin / NSeg;

  epi_raw_corr = zeros(NCol, NLin, NCha);
  for seg = 1:NSeg,

    % storage of phase correction navigator lines in phasecor1d array is
    % sparse. lines are stored consecutively along LIN dimension but are
    % assigned to different segments.
    lin0 = lines_per_segment * (seg-1) + 1;
    lin1 = lin0 + lines_per_segment - 1;

    % extract phase correction lines for this segment
    phascor1d_segment = meas.data_phascor1d(:,lin0:lin1,:,1,1,1,1,seg);

    linear_fit_coeff = mrir_artifact_ghost_compute(phascor1d_segment);

    epi_segment = epi_raw_roft(:,seg:NSeg:NLin,:,:,:,:,:,seg);
    epi_raw_corr(:,seg:NSeg:NLin,:) = mrir_artifact_ghost_correct(epi_segment, linear_fit_coeff);

  end;

  epi_img_peft = mrir_iDFT_phasencode(epi_raw_corr);
  epi_img_grid = mrir_regrid_trapezoid(epi_img_peft, meas.prot);

  epi_img_crop = mrir_image_crop(epi_img_grid, meas.prot);

  epi_rss = mrir_array_combine_rss(epi_img_crop);


  %==--------------------------------------------------------------------==%

  % reorder channel indices for phases by moving last phase to beginning of
  % list
  phases_reorder = channel_phases([end,1:end-1]);

  mode = 1;

  epi_birdcage_chan = zeros(size(epi_img_crop));

  for cha = 1:NCha,
    epi_birdcage_chan(:,:,cha) = epi_img_crop(:,:,cha) * exp( i*phases_reorder(cha) * mode);
  end;

  % sum up all birdcage channels for 1st birdcage mode
  epi_birdcage_01 = sum(epi_birdcage_chan, 3);


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_epi_birdcage_recon_kawin.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

