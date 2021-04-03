function [sigma, sigma_mag] = mrir_noise_sigma(noise);
%

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % gaussian noise model states that variance in real and imaginary part is
  % identical, so here take the average to pool and increase estimation
  % accuracy.
  sigma = [ std(real(noise(:))) + std(imag(noise(:))) ] / 2;

  if ( 0 ),
    % alternative methods:
    
    %   geometric mean
    sigma = sqrt( std(real(noise(:)))*std(imag(noise(:))) )
    
    %   arithmetic mean
    sigma = ( std(real(noise(:)))+std(imag(noise(:))) ) / 2
    
    %   pool x and y, assuming sigma is the same for each
    sigma = std([real(noise(:));imag(noise(:))])
    sigma = std(complex(real(noise(:)),imag(noise(:))))/sqrt(2)
  end;
  
  % the standard deviation measured from the magnitude data is 0.655*sigma,
  % since magnitude of noise data is rayleigh-distributed
  sigma_mag = std(abs(noise(:)));

  % estimate the percentage error based on how the measured ratio deviates
  % from theory (Haacke et al., 1999, pp. 876--877)
  scale = sqrt(2-pi/2);  %%  == 0.655
  sigma_error = ( (sigma_mag / sigma) - scale ) / scale;
  disp(sprintf('sigma_error: %2.1f%%', sigma_error*100));

  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
