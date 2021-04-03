function gfactor = mrir_array_Gfactor_1d(raw, noise, R)
%MRIR_ARRAY_GFACTOR_1D
%
% gfactor = mrir_array_Gfactor_1d(raw, noise, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/11
% $Id: mrir_array_SENSE_gfactor_1d.m,v 1.6 2007/06/20 01:46:11 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.6 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % for computations, cast data and noise as doubles to minimize roundoff
  sensitivity = mrir_sensitivity_map(double(raw), 100);
  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  
  
  Nlin = size(raw, 2); % number of phase encoded lines

  % test if the number of phase encoded lines are an integer multiple of the
  % acceleration R
  if ( mod(Nlin, R) == 0 ),
    % since Nlin divisible by R, we can exploit a faster method for this
    % special case
    gfactor = mrir_array_Gfactor_1d__special(sensitivity, covmtx, R);
  else,

    % Nlin not divisible by R: must use slower algorithm
    disp('number of PE lines not integer multiple of R! using slower, general method...');
    gfactor = mrir_array_Gfactor_1d__general(sensitivity, covmtx, R);
  end;


  return;


  
  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: