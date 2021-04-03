function [V0, C0, R0] = mrir_array_stats_sim_randmtx(Nchan)
%MRIR_ARRAY_STATS_SIM_RANDMTX
%
% [V0, C0, R0] = mrir_array_stats_sim_randcov(Nchan)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jul/03
% $Id: mrir_array_stats_sim_randmtx.m,v 1.1 2007/07/03 16:59:47 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  % generate a random covariance matrix V (must be symmetric, positive definite)
  M = complex(rand(Nchan, Nchan), rand(Nchan, Nchan));
  V0 = sqrtm(M*M');

  % due to roundoff error, imag part of diag might not == 0, so fix
  V0(logical(eye(Nchan))) = abs(diag(V0));

  % calculate corresponding correlation coefficient matrix
  w = sqrt(1./diag(V0));
  C0 = V0.*(w*w');

  % generate a random mean to produce a random correlation matrix
  mu = complex(randn(Nchan,1), randn(Nchan,1));
  R0 = V0 + [mu * mu'];


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_stats_sim_randmtx.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

