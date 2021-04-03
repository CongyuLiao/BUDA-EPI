function sensitivity = mrir_sensitivity_map(varargin)
%MRIR_SENSITIVITY_MAP  estimate coil sensitivity maps from image data
%
% sensitivity = mrir_sensitivity_map(raw, [smooth_method], [smooth_param], [noise_std], [reference]);
%
% default: no smoothing, "smooth_method" = 0

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/28
% $Id: mrir_sensitivity_map.m,v 1.1 2007/06/01 17:44:07 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  if ( iscell(varargin{1}) ),
    raw  = varargin{1}{1};
    prot = varargin{1}{2};
  else,
    raw  = varargin{1};
    prot = read_meas_prot__struct;
  end;

  smooth_method = 0;
  if ( nargin >= 2 ),
    smooth_method = varargin{2};
  end;

  smooth_param = 1.0;
  if ( nargin >= 3 ),
    smooth_param = varargin{3};
  end;

  noise_std = [];
  if ( nargin >= 4 ),
    noise_std = varargin{4};
  end;

  FLAG__REF_NORMALIZE = 0;
  raw_reference = [];
  if ( nargin >= 5 ),
    FLAG__REF_NORMALIZE = 1;
    raw_reference = varargin{5};
  end;

  DEBUG = 0;


  %==--------------------------------------------------------------------==%

  % method1: truncate and zero-pad phase encoded lines (allowing ringing)
  % method2: truncate and apply Hann window apodization symmetrically
  % method3: smooth with symmetric Gaussian kernel
  % method4: fit with an n-degree polynomial
  % method5: fit with a Savitzky-Golay filter
  % method6: fit with an exponential
  % method7: fit with b-splines

  % - masking based on threshold from measured noise variance
  % - normalization with body coil
  % - normalization with sum-of-squares

  Nfreq = size(raw, 1);
  Nphas = size(raw, 2);
  Nchan = size(raw, 3);
  Npart = size(raw, 9);
  Nslic = size(raw,10);

  fft_normalize = Nfreq * Nphas * Npart;


  %==--------------------------------------------------------------------==%
  %%% reconstruct sensitivity images and reference images

  % conventional reconstruction (asssuming no EPI sensitivity maps!)
  sensitivity_img = mrir_conventional(double(raw), prot);


  if ( FLAG__REF_NORMALIZE ),

    if ( isempty(raw_reference) ),
      % use sum-of-squares combo for reference (magnitude only!)
      reference = mrir_array_combine_rss(sensitivity_img);
    end;

    reference = mrir_conventional(raw_reference, prot);

    if ( size(raw_reference, 3) > 1 ),
      disp('combining multi-channel reference data for magnitude-only normalization');
      reference = mrir_array_combine_rss(reference);
    end;

  else,

    reference = 1;

  end;


  %==--------------------------------------------------------------------==%
  %%% smooth sensitivity (and reference) images

  switch lower(smooth_method),
    %-------------------------------------------------------------------*
   case {0, 'none'},
    %    disp('skipping sensitivity map spatial smoothing...');
    
    sensitivity_smoothed = sensitivity_img;

    
    %-------------------------------------------------------------------*
   case {1, 'pad'},
    percentage = smooth_param;

    % (safely) assume that resolution is limited in the phase encoding direction
    extent_from_center = round((Nphas * percentage)/2);
    zero_ind_lo = 1:(Nphas/2 - extent_from_center);
    zero_ind_hi = (Nphas/2 + extent_from_center + 1):Nphas;

    raw_zeropadded = double(raw);

    raw_zeropadded(:,zero_ind_hi,:) = 0;
    raw_zeropadded(:,zero_ind_lo,:) = 0;

    % assume no ramp sampling
    sensitivity_smoothed = mrir_conventional(raw_zeropadded, prot);


    %-------------------------------------------------------------------*
   case {2, 'hann'},

    percentage = smooth_param;
    sensitivity_smoothed = mrir_image_smooth__kspace_subset(sensitivity_img, percentage);


    %-------------------------------------------------------------------*
   case {3, 'gauss'},

    percentage = smooth_param;

    apod_FWHM_pixel = 2* 1/percentage;
    apod_std_pixel = apod_FWHM_pixel / (2*sqrt(2*log(2)));

    t_PE = -(Nphas/2 - 1):(Nphas/2);
    apod_PE = exp(-0.5*(t_PE.^2) / (apod_std_pixel^2));
    %  apod_PE = sinc(t_PE/2/pi).*hamming(Nphas).';
    Apod_PE = abs(ifftshift(fft(fftshift(apod_PE / sum(apod_PE))))'');
    APOD_PE = repmat(Apod_PE, [Nfreq, 1]);

    t_FE = -(Nfreq/2 - 1):(Nfreq/2);
    apod_FE = exp(-0.5*(t_FE.^2) / (apod_std_pixel^2));
    %  apod_FE = sinc(t_FE/2/pi).*hamming(Nfreq).';
    Apod_FE = abs(ifftshift(fft(fftshift(apod_FE / sum(apod_FE)))).');
    APOD_FE = repmat(Apod_FE, [1, Nphas]);

    Apod = Apod_FE*Apod_PE * sqrt(2);
    %  Apod = Apod / sum(Apod(:));

    for ind = 1:Nchan,
      raw_smoothed(:,:,ind) = double(raw(:,:,ind)).*Apod;
    end;

    sensitivity_smoothed = mrir_conventional(raw_smoothed);


    if ( DEBUG ),
      window = abs(t_PE-0.5) <= Nphas/2*percentage;

      figure('name', mfilename); plot(t_PE, abs(Apod), 'b'); hold on
      plot(t_PE, window, 'r'); axis tight
    end;


    %-------------------------------------------------------------------*
   case {4, 'poly'},
    %-------------------------------------------------------------------*
   case {5, 'golay'},
    %-------------------------------------------------------------------*
   case {6, 'exp'},
    %-------------------------------------------------------------------*
   case {7, 'spline'},

    %-------------------------------------------------------------------*
   otherwise,
    disp('method "%s" not recognized', num2str(smooth_method));
  end;





  %==--------------------------------------------------------------------==%
  %%% normalized smoothed sensitivity maps with smoothed reference images

  % normalize with reference data
  for chan = 1:Nchan,
    sensitivity(:,:,chan) = sensitivity_smoothed(:,:,chan) ./ reference;
  end;

%%  if ( isempty(prot.flReadoutOSFactor) ),
%%    sensitivity = sensitivity / sqrt(2);
%%  else,
%%    sensitivity = sensitivity / sqrt(prot.flReadoutOSFactor);
%%  end;

  %==--------------------------------------------------------------------==%


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sensitivity_map.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:






  %==--------------------------------------------------------------------==%
  %%% generate signal/background mask from (unsmoothed) sensitivity map

  img_comb = mrir_array_combine_rss(sensitivity_img) / fft_normalize;


  % since only have one measurement, compute standard deviation of image
  % intensities over a sliding 5x5 window
  std_comb = sqrt(ordfilt2(img_comb, 2, ones(5,5)));

  %avg_comb = filter2(ones(5,5)/25, img_comb, 'same');
  %avg_norm = (avg_comb - min(avg_comb(:))) / (max(avg_comb(:)) - min(avg_comb(:)));

  %sigma_data = std_comb*scale_factor/sqrt(2-pi/2);
  %sigma_data = avg_comb*scale_factor/sqrt(pi/2);

  if ( isempty(noise_std) ),

    threshold_percent = 20;
    threshold_value = mean(mean(img_comb)) * threshold_percent/100;
    mask_raw = img_comb > threshold_value;

  else,
    noise_mag_std = noise_std * sqrt(2-pi/2);

    % signal is above 4 noise standard deviations
    threshold_std = 4;
    threshold_value = noise_mag_std * threshold_std;
    mask_raw = std_comb > threshold_value;

  end;

  mask_padded = padarray(mask_raw, [Nfreq/2, Nphas/2]);
  step1 = bwmorph(mask_padded, 'clean', inf);
  step2 = bwmorph(step1, 'fill',        inf);
  step3 = bwmorph(step2, 'dilate',      min([Nfreq,Nphas])/2-2);
  step4 = bwmorph(step3, 'erode',       min([Nfreq,Nphas])/2-2);

  mask = step4;






  kernel_boundary = strel('disk', 1);

  % strip off one pixel around boundary since usually find artifact there
  mask_trim = imerode(step4, kernel_boundary);

  % extract boundary pixels
  mask_interior = imerode(mask_trim, kernel_boundary);
  mask_boundary = (mask_trim - mask_interior) == 1;

  kernel_extend = strel('disk',   ceil(3*apod_std_pixel));
  kernel_smooth = strel('disk', 3*ceil(3*apod_std_pixel));
  kernel_square = strel('square', 1);

  mask_extend = (imdilate(mask_trim, kernel_extend) - mask_trim) == 1;


  figure('name', mfilename);
  imagesc( 2*mask_trim + mask_extend ); axis image; colormap(jet(5)); colorbar;

  figure('name', mfilename);
  imagesc( 2*mask_trim - mask_boundary ); axis image; colormap(jet(5)); colorbar;



  %     [ind_boundary_row, ind_boundary_col] = find(mask_boundary);
  %     [ind_extend_row, ind_extend_col] = find(mask_extend);
  %     [ind_trim_row, ind_trim_col] = find(mask_trim);
  %
  %     for rr = 1:length(ind_boundary_row),
  %       for cc = 1:length(ind_boundary_col),
  %
  %         dist = sqrt( (ind_trim_row - ind_boundary_row(rr)).^2 ...
  %                    + (ind_trim_col - ind_boundary_col(cc)).^2 );
  %
  %
  %       end;
  %     end;
  %
  %     figure('name', mfilename);
  %     imagesc(img_smooth); axis image;
  %
  %
  %       img_smooth = imdilate(img_boundary, kernel_boundary);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_square);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_boundary);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_square);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_boundary);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_square);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_boundary);
  % figure; imagesc(img_smooth);
  %       img_smooth = imdilate(img_smooth, kernel_boundary);
  %
  %


  h = fspecial('gaussian', ceil(3*apod_std_pixel), apod_std_pixel);

  img_output = zeros(size(img_crop));
  for ind = 1:Nchan,

    img_chan = abs(img_crop(:,:,ind));

    img_boundary = zeros(size(img_chan));
    img_boundary(mask_boundary) = img_chan(mask_boundary);

    img_smooth = imdilate(img_boundary, kernel_smooth);

    img_extend = zeros(size(img_chan));
    img_extend(mask_trim) = img_chan(mask_trim);
    %      img_extend(mask_extend) = img_smooth(mask_extend);
    %      img_extend(mask_extend) = mean(img_chan(mask_trim));
    img_max = maxfilt2(img_boundary, ceil(3*apod_std_pixel)*[1,1]+1);
    img_extend(mask_extend) = img_max(mask_extend);

    img_conv = conv2(img_extend, h, 'same') .* exp(i*angle(img_crop(:,:,ind)));

    img_mask = zeros(size(img_conv));
    img_mask(mask) = img_conv(mask);

    img_output(:,:,ind) = img_mask;

  end;

