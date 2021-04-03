function i_full = mrir_partial_fourier(h_part, prot, varargin)
%MRIR_PARTIAL_FOURIER  partial fourier reconstruction in phase encoding direction
%
% h_full = mrir_partial_fourier(h_part, prot)
% h_full = mrir_partial_fourier(h_part, prot, METHOD)
% h_full = mrir_partial_fourier(h_part, prot, METHOD, KSPACE_FILL)
%
% note that if k-space is partially covered in more than one dimension, the
% homodyne processing method can only be used for one dimension, and must
% follow zero-filling reconstructions along all other partially covered
% dimensions.
%
% fourier transform must always occur after zero-filling in any dimension
% before zero-filling (or homodyne) can be applied in other dimensions, so
% to force this behavior this function also computes fourier transform as
% the last step.
%
% rules summary:
%  1) all fully-acquired dimensions must be iDFT'd before any PF processing
%  2) if multiple PF dimensions, only last to be processed can use homodyne
  
% references:
%
%  Bernstein et al, section 13.4.2
%
%  Pauly course notes, 2003, chapter 2.

% examples:
%
%   data = getfield(load('mri', 'D'), 'D');
%   i_orig = double(squeeze(data(:,:,1,20)));
%   prot.iNoOfFourierColumns = 128;
%   prot.lPhaseEncodingLines = 128;
%   h_orig = mrir_fDFT_phasencode(i_orig);
%   i_test1 = mrir_partial_fourier(h_orig(:,32+1:end), prot, 'homodyne', 'pre');
%   i_test2 = mrir_partial_fourier(h_orig(:,1:end-32), prot, 'homodyne', 'post');
%
%   k_orig = mrir_fDFT_freqencode(h_orig);
%   k_part = k_orig(32+1:end,32+1:end);
%   h_test = mrir_partial_echo(k_part, prot, 'zero-fill');
%   h_test = mrir_partial_fourier(k_part.', prot, 'zero-fill').';
%   i_test3 = mrir_partial_fourier(h_test, prot, 'zero-fill');
%   i_test4 = mrir_partial_fourier(h_test, prot, 'homodyne');
%
%   h_fail = mrir_partial_fourier(k_part.', prot, 'homodyne').';
%   i_fail1 = mrir_partial_fourier(h_fail, prot, 'zero-fill');
%   i_fail2 = mrir_partial_fourier(h_fail, prot, 'homodyne');
  
% TODO: generalize so that handles partial fourier in frequency- and
% partition-encoding directions for both methods

% TODO: adapt "read_meas_dat" to determine automatically whether fill should
% be on the positive or negative side of the k-space origin (i.e., whether
% "KSPACE_FILL" is 'pre' or 'post').
  
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/16
% $Id: mrir_partial_fourier.m,v 1.2 2007/10/28 04:29:15 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global FLAG__FIRST_CALL; if ( isempty(FLAG__FIRST_CALL) ), FLAG__FIRST_CALL = 0; end;

  
  % default method is zero filling
  METHOD = 'zero-fill';
  if ( nargin >= 3 ),

    % options are 'zero-fill' or 'homodyne'
    METHOD = varargin{1};
  
  end;


  % default omitted k-space data occurs BEFORE echo so TE is shorter in EPI,
  % but sometimes is after echo to shorten total acquisition time in
  % conventional fourier imaging.
  KSPACE_FILL = 'pre';
  if ( nargin >= 4 ),

    % options are 'pre' or 'post'
    KSPACE_FILL = varargin{2};
  
    if ( ~ismember(KSPACE_FILL, {'pre', 'post'}) ),
      error('unrecognized parameter value "%s"; fill region must be either "pre" or "post"', KSPACE_FILL);
    end;
  
  end;
  
  
  %==--------------------------------------------------------------------==%

  Nlin_full = prot.lPhaseEncodingLines;  % also, "Nlin_part / prot.ucPhasePartialFourier"

  if ( Nlin_full == mrir_ice_dimensions(h_part, 'lin') ),

   % disp(sprintf('==> [%s]: no partial fourier in LIN dimension detected, computing iDFT directly...', mfilename));
    i_full = mrir_iDFT_phasencode(h_part);

    return;
  end;


  %==--------------------------------------------------------------------==%

  switch lower(METHOD(1)),
   case 'z',   % zero-fill
    i_full = mrir_partial_fourier__zerofill(h_part, Nlin_full, KSPACE_FILL);
   
    % compensate intensity scale for true number of samples
    i_full = i_full / sqrt(prot.ucPhasePartialFourier);
  
   case 'h',   % homodyne,

    if ( FLAG__FIRST_CALL ),
      disp(sprintf('==> [%s]: reconstructing using homodyne processing...', mfilename));
    end;
      
    i_full = mrir_partial_fourier__homodyne(h_part, Nlin_full, KSPACE_FILL);
   otherwise,
    error('unknown partial fourier recon method: "%s"', METHOD);
  end;


  return;



%**************************************************************************%
function i_full = mrir_partial_fourier__zerofill(h_part, Nlin_full, KSPACE_FILL)


  Nlin_part = mrir_ice_dimensions(h_part, 'lin');

  pad_lines = Nlin_full - Nlin_part;


  %                          1         2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  %h_full = padarray(h_part, [0 pad_lines 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 0, KSPACE_FILL);
  h_full = mrir_zeropad( h_part,[0 pad_lines 0 0 0 0 0 0 0 0 0 0 0 0 0 0],KSPACE_FILL);

  i_full = mrir_iDFT_phasencode(h_full);


  return;


%**************************************************************************%
function i_full = mrir_partial_fourier__homodyne(h_part, Nlin_full, KSPACE_FILL)

  dims_part = size(h_part);

  Nlin_part = mrir_ice_dimensions(h_part, 'lin');


  index_cent = Nlin_full / 2;
  Nlin_symm = (Nlin_part - index_cent) * 2;

  index_symm = (Nlin_full - Nlin_part + 1) : (Nlin_full - Nlin_part + 1) + Nlin_symm - 1;

  index_zero = index_symm(1) - 1;

  % ramp part of pre-weighting function. calculate so that corners of
  % weighting occur directly on data outside of symmetric region, i.e., none
  % of the symmetric region will have a zero weight.
  preweight_ramp = linspace(0, 2, Nlin_symm+2);
  preweight_ramp = preweight_ramp(2:end-1);

  preweight = [zeros(index_zero,1); preweight_ramp'; 2*ones(index_zero,1)].';

  
  switch KSPACE_FILL,
   case 'pre',
   case 'post',
    preweight = fliplr(preweight);
  end;

  
  dims_full = dims_part;
  dims_full(2) = Nlin_full;

  h_symm = zeros(dims_full);

  switch KSPACE_FILL,
   case 'pre',
    %      1,         2,3,4,5,6,7,8,9,0,1,2,3,4,5,6           1,                      2,3,4,5,6,7,8,9,0,1,2,3,4,5,6
    h_symm(:,index_symm,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = h_part(:,index_symm - index_zero,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
   case 'post',
    %      1,         2,3,4,5,6,7,8,9,0,1,2,3,4,5,6           1,                      2,3,4,5,6,7,8,9,0,1,2,3,4,5,6
    h_symm(:,index_symm,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = h_part(:,index_symm,             :,:,:,:,:,:,:,:,:,:,:,:,:,:);
  end;
  
  % image of symmetric part
  i_symm = mrir_iDFT_phasencode(h_symm, 'lin');


  % hybrid space of anti-symmetric part
  h_anti = zeros(dims_full);

  switch KSPACE_FILL,
   case 'pre',
    %      1,                2,3,4,5,6,7,8,9,0,1,2,3,4,5,6           1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6
    h_anti(:,index_symm(1):end,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = h_part(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
   case 'post',
    %      1,                2,3,4,5,6,7,8,9,0,1,2,3,4,5,6           1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6
    h_anti(:,1:index_symm(end),:,:,:,:,:,:,:,:,:,:,:,:,:,:) = h_part(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
  end;
  

  % high-pass filtered version in hybrid space
  h_high = h_anti .* repmat(preweight, [dims_full(1), 1, dims_full(3:end)]);

  % high-pass filtered version in image space
  i_high = mrir_iDFT_phasencode(h_high, 'lin');

  % rotate phases in high-pass filtered image with symmetric part image to reach final result.
  i_full = real( i_high .* exp(-i*angle(i_symm)) );


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_partial_fourier.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:




%%%  % apodized high-pass filter (bernstein2004:handbook, pg. 551)
%%%  wid = 2;
%%%  homodyne_highpass2 = [repmat(0,start-w/2,1); repmat(1,symmetric+2,1); repmat(2,start-1,1)];
%%%
%%%  cos(pi*(abs(k) - (k0 - w/2))/2/w).^2 + 0;
%%%  cos(pi*(abs(k) - (k0 + w/2))/2/w).^2 + 1;
