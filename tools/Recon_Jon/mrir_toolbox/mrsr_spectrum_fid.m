function [fid_spectrum, f_Hz] = mrir_spectrum_fid(fid, dwelltime_us)
%

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/19
% $Id: mrir_spectrum_fid.m,v 1.1 2007/09/23 00:30:24 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  DO__PLOT = 1;


  fid_avg = mean(fid, 4);

  fid_spectrum = fft(fid_avg, [], 1);


  fid_amplitude_max = max(abs(fid_spectrum))



  delta_t = dwelltime_us / 1e6 / mrir_ice_dimensions(fid, 'COL');  % sec
  f_max = 1/delta_t;    % Hz

  w = linspace(-pi, +pi, mrir_ice_dimensions(fid_spectrum, 'COL')+1);

  f = w / 2 / pi * f_max;
  f(end) = [];

  if ( DO__PLOT ),

    figure('name', mfilename);
    plot(f, fftshift(abs(squeeze(fid_spectrum))))
    set(gca, 'XDir', 'reverse');
    xlabel('frequency (Hz)');
    ylabel('FID amplitude spectrum');

    if ( ~isempty(inputname(1)) ),
      set(title(inputname(1)), 'Interpreter', 'none');
    end;

  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_spectrum_fid.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
