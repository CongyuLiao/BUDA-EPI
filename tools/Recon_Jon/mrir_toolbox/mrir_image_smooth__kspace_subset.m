function sensitivity_smoothed = mrir_image_smooth__kspace_subset(sensitivity_img, percentage)

% mrir_image_smooth__kspace_subset(img, percentage)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/30
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  
  Nfreq = size(sensitivity_img, 1);
  Nphas = size(sensitivity_img, 2);
  Nchan = size(sensitivity_img, 3);
  Npart = size(sensitivity_img, 9);
  Nslic = size(sensitivity_img,10);

  fft_normalize = Nfreq * Nphas * Npart;

  

    if ( percentage > 1 ),
      % user did not pass percentage correctly!
      percentage = percentage / 100;
    end;
    
    % force filter to be an integer and a multiple of 2
    Nfilter = ceil( Nphas*percentage )*2;

    % generate 2D Hann window
    Hann1 = hanning(Nfilter);
    Hann2 = repmat(Hann1, [1, Nfilter]) .* repmat(Hann1', [Nfilter, 1]);

    % normalize for unity noise gain from filtering (Kellman & McVeigh, 2005)
    Hann2 = Hann2 / sqrt(mean(mean(Hann2.^2)));

    % zero-pad to twice the image size (to avoid wrapping after filtering)
    H = padarray(Hann2, (2*Nphas-Nfilter)/2 * [1 1]);
    sensitivity_padded = padarray(sensitivity_img / fft_normalize / 2, Nphas/2 * [1 1]);

    % transform sensitivity map back into k-space for windowing (to avoid
    % heavy burden of convolution!)
    S = mrir_fDFT_phasencode(mrir_fDFT_freqencode(sensitivity_padded));

    % window channel-by-channel
    for chan = 1:Nchan,
      S_apod(:,:,chan) = S(:,:,chan) .* H;
    end;

    % transform back into image space
    s_apod = mrir_iDFT_phasencode(mrir_iDFT_freqencode(S_apod));

    % extract out center from padded image
    ind_unpad = [1:Nphas] + (Nphas/2);
    sensitivity_smoothed = s_apod(ind_unpad, ind_unpad, :);

    %1/sqrt(percentage);

    scaling = max(max(abs(sensitivity_smoothed(:,:,1))));
    
    if ( scaling ~= 0 ),
      sensitivity_smoothed = sensitivity_smoothed * max(max(abs(sensitivity_img(:,:,1)))) / max(max(abs(sensitivity_smoothed(:,:,1))));
    end;
    
    
    return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
