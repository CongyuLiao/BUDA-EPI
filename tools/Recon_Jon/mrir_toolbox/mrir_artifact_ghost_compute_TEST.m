function linear_fit_coeff = mrir_artifact_ghost_compute(phascor)
%MRIR_ARTIFACT_GHOST_COMPUTE
%
% phase_correct_odd = mrir_artifact_ghost_compute(phascor_raw)
%
% returns a 2 x Nphascor matrix of linear fit coefficients, where Nphascor
% is the number of unique phase correction lines collected at the center of
% k-space (e.g., the number of coil channels). at least one pair of phase
% correction navigators must be collected, and hopefully this measurement is
% repeated several times. the first row is the offset term, and the second
% row is the slope. all units are radians.

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2006/dec/31
% $Id: mrir_artifact_ghost_compute.m,v 1.2 2007/11/14 03:05:15 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG;   if ( isempty(DEBUG)   ), DEBUG   = 0; end;
  global VERBOSE; if ( isempty(VERBOSE) ), VERBOSE = 0; end;

  
  %--------------------------------------------------------------------------%

  % ALGORITHM OUTLINE:
  %
  % (1) average (complex-valued) subset of odd and even lines
  % (2) subtract phases of even line from odd line FIRST
  % (3) fit a first-order polynomial to the phase difference

  % perform IDFT along readout direction
  projections = mrir_iDFT_freqencode(phascor);

  % crop out extra samples due to oversampling k-space
  projections = mrir_image_crop(projections);
  projections = mrir_image_crop(projections);

  dims = size(projections);

  % pick the first several lines (before too much T2* delay)
  %                                1 2 3 4 5 6 7       8 9 0 1 2 3 4 5 6
  Nref1 = min([8, size(projections(:,:,:,:,:,:,:,1:2:end,:,:,:,:,:,:,:,:),2)]);
  Nref2 = min([8, size(projections(:,:,:,:,:,:,:,2:2:end,:,:,:,:,:,:,:,:),2)]);

  % separate out odd and even reference lines (assuming direction
  % reverses every line)
  % NOTE: convention is that first refline is a reflect (odd) and second is not (even)!
  %                        1       2 3 4 5 6 7       8 9 0 1 2 3 4 5 6
  ref_o = mean(projections(:,1:Nref1,:,:,:,:,:,1:2:end,:,:,:,:,:,:,:,:), 2);  % line "0" is odd
  ref_e = mean(projections(:,1:Nref2,:,:,:,:,:,2:2:end,:,:,:,:,:,:,:,:), 2);  % line "1" is even

  dims(8) = mrir_ice_dimensions(projections, 'seg') / 2;
  
  
  % collapse even and odd reference lines into a Nreadout x Nlines array
  ref_e = reshape(ref_e, [dims(1), prod(dims(3:end))]);
  ref_o = reshape(ref_o, [dims(1), prod(dims(3:end))]);

  % [Nreadout x Nlines]
  phasediff = unwrap( angle(ref_e(:,:)./ref_o(:,:)), [], 1);

  % simple test for 2pi shifts
  % [       1 x Nlines]
  phasediff_center_pos = (mean(phasediff, 1) > +pi) * 2*pi;
  % [Nreadout x Nlines]
  phasediff_offset_pos = repmat(phasediff_center_pos, size(phasediff,1), 1);

  phasediff_center_neg = (mean(phasediff, 1) < -pi) * 2*pi;
  phasediff_offset_neg = repmat(phasediff_center_neg, size(phasediff,1), 1);

  phasediff = phasediff - phasediff_offset_pos;
  phasediff = phasediff + phasediff_offset_neg;

  % keep only the center 50% of samples, since noise at ends can bias fit
  phasediff = phasediff(ceil(1*end/4):floor(3*end/4),:);
  
  
  % build indices into readout direction for offset and slope of phasediff
  % (since "phasediff" is truncated, these will not be same size as the data
  % correction terms)
  x_refline = [ones(size(phasediff,1),1), ([0:size(phasediff,1)-1]-size(phasediff,1)/2)'];

  % compute linear (1st order poly) fit for each reference line
  linear_fit_coeff = pinv(x_refline)*phasediff;

  linear_fit_coeff = [linear_fit_coeff(1,:); linear_fit_coeff(2,:)];


  coeff_avg = mean(linear_fit_coeff, 2);

  if ( VERBOSE || DEBUG ),
    disp(sprintf('<i> [%s]: average parameters:  offset = %2.2f deg; slope = %2.2f deg', ...
		 mfilename, rad2deg(coeff_avg)));
  end;
  
  if ( DEBUG ),
    disp(sprintf('<d> [%s]: generating debugging figures', mfilename));
    mrir_artifact_ghost_compute__DEBUG(phascor, phasediff, linear_fit_coeff, ref_e, ref_o, dims);
  end;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_artifact_ghost_compute.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
