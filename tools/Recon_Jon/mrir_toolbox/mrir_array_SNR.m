function snr = mrir_array_SNR(raw, noise, METHOD)
%MRIR_ARRAY_SNR
%
% snr = mrir_array_SNR(raw, noise, METHOD)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( ~ismember(METHOD, {'opt', 'rss', 'ind'} ) ),
    error('unknown SNR calculation method: "%s"', METHOD);
  end;
  
  
  %==--------------------------------------------------------------------==%

  % for computations, cast data and noise as doubles to minimize roundoff
  sensitivity = mrir_sensitivity_map(double(raw), 100);
  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  
  switch(lower(METHOD)),
   case 'opt',
    snr = mrir_array_SNR_opt(raw, noise);
   case 'sos',
    snr = mrir_array_SNR_rss(raw, noise);
   case 'ind',
    snr = mrir_array_SNR_ind(raw, noise);
   otherwise,
    error('unsupported SNR calculation method: "%s"', METHOD);
  end;
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
