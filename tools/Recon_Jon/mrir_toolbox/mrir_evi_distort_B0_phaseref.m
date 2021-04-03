function evi_undistort = mrir_evi_distort_B0_phaseref(evi_pelinft, phaseref_raw)
%MRIR_EVI_DISTORT_B0_PHASEREF
%
% evi_undistort = mrir_evi_distort_B0_phaseref(evi_pelinft, phaseref_raw)

% van der Zwaag W, Francis S, Bowtell R. "Improved echo volumar imaging
% (EVI) for functional MRI". Magn Reson Med. 2006 Dec;56(6):1320-7. PMID:
% 17089364.


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/24
% $Id: mrir_evi_distort_B0_phaseref.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  phaseref_raw_roft = mrir_iDFT_freqencode(phaseref_raw);
  phaseref_img_peft = mrir_iDFT_phasencode(phaseref_raw_roft, 'lin');

  % phase map extracted from field map after Fourier transformation with
  % respect to frequency-encoded and primary phase-encoded directions.
  phasemap = angle(phaseref_img_peft);


  % remove phase
  evi_undistort = evi_pelinft .* exp(-i*phasemap);

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_distort_B0_phaseref.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
