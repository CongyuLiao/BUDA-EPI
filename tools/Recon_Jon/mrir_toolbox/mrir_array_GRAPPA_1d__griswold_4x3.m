function k_raw = mrir_array_GRAPPA_1d__griswold_4x3(k_dat, k_acs, evp)
%mrir_array_GRAPPA_1d__griswold_4x3  wrapper around "opengrappa_4x3"
%
% k_raw = mrir_array_GRAPPA_1d__griswold_4x3(k_dat, k_acs, evp)

% example:
%
%   [k_dat, k_acs] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);
%
%   k_raw = mrir_array_GRAPPA_1d__griswold_4x3(k_dat, k_acs, meas.evp)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/02
% $Id: mrir_array_GRAPPA_1d__griswold_4x3.m,v 1.1 2008/11/03 03:08:38 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  sig = permute(k_dat, [3, 2, 1]);
  acs = permute(k_acs, [3, 2, 1]);

  af = evp.NAFLin;


  [recon, sigrecon, ws] = opengrappa_4x3(sig, acs, af);


  k_raw = permute(sigrecon, [3, 2, 1]);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_1d__griswold_4x3.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: