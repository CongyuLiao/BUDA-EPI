function rho = mrir_array_SENSE_gfactor_global_corrcoef(g)
%MRIR_ARRAY_SENSE_GFACTOR_GLOBAL_CORRCOEF  global correlation coefficient map
%
% rho = mrir_array_SENSE_gfactor_global_corrcoef(g)
%
%
% James, F. (2006). Statistical methods in experimental physics. Hackensack,
% NJ: World Scientific, pg. 28

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/03
% $Id: mrir_array_SENSE_gfactor_global_corrcoef.m,v 1.1 2008/01/03 23:18:32 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  rho = sqrt(1 - 1./(g.^2));

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_global_corrcoef.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
