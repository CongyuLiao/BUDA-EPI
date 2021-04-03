function snr = mrir_array_SNR_timeseries(img_timeseries)
%MRIR_ARRAY_SNR_TIMESERIES
%
% snr = mrir_array_SNR_timeseries(img_timeseries);

% convention: assume repetitions is in either dimension 7 (following ICE),
% or, if fewer than 7 dimensions, in the last dimension.


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/20
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dim = min([ ndims(img_timeseries), 7 ]);

  % there is an asymmetry in how MATLAB computes mean and std of
  % complex-valued data. the mean will be a complex number whose real part
  % is the mean of the real part of the data, and whose imaginary part is
  % the mean of the imaginary part of the data. however, the standard
  % deviation is always a real number corresponding to the standard
  % deviation of the magnitude of the data!


  % cases to consider:
  %  1) input is the abs of complex-valued data, e.g., rss or opt
  %  2) input is real-valued and non-negative
  %  3) input is real-valued containing positive and negative values

  % normalize by sqrt(2) to compensate for real and imagingary parts
  % combined in the magnitude data:
  %%% mean(z) == m + i*m, abs(mean(z)) == sqrt(m^2 + m^2) == m * sqrt(2)

  % NOTE: this model is FALSE. typically signal *is* the image magnitude,
  % and the real and imaginary parts are not the same value. there should
  % therefore never be any sqrt(2) normalization. however, whenever
  % taking the magnitude 
  
  if ( isreal(img_timeseries) ),
    NORMALIZE_MAG = sqrt(2);
  else,
    NORMALIZE_MAG = 1;
  end;
  
  
  img_timeseries_avg = mean( img_timeseries, dim) / NORMALIZE_MAG;
  img_timeseries_std =  std( img_timeseries, 0, dim );


  % note that if "img_timeseries" is actually real-valued, then the
  % arithmetic average of "avg" and "std" with both be off by a factor of
  % two, which will cancel in the snr ratio!
  snr = abs(img_timeseries_avg) ./ img_timeseries_std;

  % complex-valued SNR does not make sense.
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
