function [dphi0, dphi1m, dphi1, error_rms] = mrir_artifact_ghost_compute__PC_algo_siemens_local(phascor)
%MRIR_ARTIFACT_GHOST_COMPUTE__PC_ALGO_SIEMENS_local  following algorithm in AdvOnlineTse
%
% [phi0_seg1, phi1_seg1, phi1_seg2] = mrir_artifact_ghost_compute__PC_algo_siemens_local(phascor)
%
%
% see also MRIR_ARTIFACT_GHOST_COMPUTE__PC_ALGO_SIEMENS_VB.

% based on the methods described in:
%
% Feiweier T. Magnetic Resonance Method and Apparatus to Determine Phase
% Correction Parameters. U.S. Patent No. 8,497,681, Issued 7/30/2013.

% there are four essential elements of this method:
%
% 1) the "local" calculation in which the dot product is weighted by the
% signal amplitudes
%
% 2) the additional "weighting function", selectable by the user, to specify
% which sub-regions of the image are to be prioritized in the phase
% correction calculation, e.g., a hanning window centered on some anatomy
%
% 3) the allowance for higher-order polynomials or "power series"
%
% 4) the sliding-window approach or segmented approach to breaking up the
% data into smaller chunks, each with their own fitting

% as of the first version of this function, only element #1 has been
% implemented, because it is the most straightforward (i.e., it doesn't
% require any parameters or decisions) and because it appears to do most of
% the work---which is presumably why siemens has named the new method "local
% phase correction"

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/23
% $Id: mrir_artifact_ghost_compute__PC_algo_siemens_local.m,v 1.1 2015/08/23 19:23:04 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % perform IDFT along readout direction
  hyb = mrir_iDFT_freqencode(phascor);

  %%% cropping should not be needed in the local phase correction since the
  %%% weighting is built in
  % crop out extra samples due to oversampling k-space
  %  hyb = mrir_image_crop(hyb);
  %  hyb = mrir_image_crop(hyb);

  Ncol = size(hyb, 1);
  Ncha = size(hyb, 3);


  % first-order correction (auto-cor), then zeroth-order correction (cross-cor)


  % notation following Feiweier (2013)

  % NOTE: convention is that segment "0" is a normal line and segment "1" is
  % a reflected line -- first or second line varies, but reflected always
  % segments "1 + 2n"

  % the first phase correction data set is determined by averaging the 1st
  % and 3rd lines, as in Heid (2000)

  ind1 = find(sum(sum(sum(sum(hyb(:,:,:, 1,1,1,1, 01, :,:), 1), 3), 9), 10));
  ind2 = find(sum(sum(sum(sum(hyb(:,:,:, 1,1,1,1, 02, :,:), 1), 3), 9), 10)); 

  %              1  2   3 4 5 6 7 8 9 0 1 2 3 4 5 6
%   S1  = mean(hyb(:,ind2,:,:,:,:,:,2,:,:,:,:,:,:,:,:), 2);  %% NOTE this is flipped for now!
%   S2  = mean(hyb(:,ind1,:,:,:,:,:,1,:,:,:,:,:,:,:,:), 2);
   S1  = mean(hyb(:,ind2,:,:,:,:,:,1,:,:,:,:,:,:,:,:), 2);  %% unflipped
  S2  = mean(hyb(:,ind1,:,:,:,:,:,2,:,:,:,:,:,:,:,:), 2);

  n = 1:(Ncol-1);

  % these are shifted versions of one another (but calculate both to follow
  % notation of Feiweier)
  C_xnp0 = S2(n+0, 1, :, :,:,:,:,:,:,:) .* conj(S1(n+0, 1, :, :,:,:,:,:,:,:));
  C_xnp1 = S2(n+1, 1, :, :,:,:,:,:,:,:) .* conj(S1(n+1, 1, :, :,:,:,:,:,:,:));

  f = 1;

  % equation [4]: fitting the difference S2./S1 or seg1 - seg2
  dphi1 = angle(sum(f .* C_xnp1 .* conj(C_xnp0), 1));

  % for backwards compatibility
  dphi1m = zeros(size(dphi1));

  % coords in hybrid space
  x = [([0:size(hyb,1)-1]-size(hyb,1)/2)'];

  dims = size(S1+S2);

  X = repmat(x, [1, dims(2:end)]);
  Dphi1 = repmat(dphi1, [dims(1), ones(1, length(dims)-1)]);

  g = 1;

  dphi0 = angle(sum( g .* S1 .* conj(S2) .* exp(i*X.*Dphi1), 1));

  %=---
  
  Dphi0 = repmat(dphi0, [dims(1), ones(1, length(dims)-1)]);

  error_rms = angle(mean(S1 .* conj(S2) .* exp(i*X.*Dphi1) .* exp(-i*Dphi0)));
  
  
  % [dphi0, ~, dphi1, error_rms] = mrir_artifact_ghost_compute__PC_algo_siemens_local(phascor)
  
  
%  for cha = 1:Ncha,
%
%    phi0(:, 1, cha) = angle(sum( g .* S1(:, 1, cha) .* conj(S2(:, 1, cha) ) .* vec(exp(i*(j - (Ncol/2))*dphi1(:,:,cha))), 1));
%
%  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_artifact_ghost_compute__PC_algo_siemens_local.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
