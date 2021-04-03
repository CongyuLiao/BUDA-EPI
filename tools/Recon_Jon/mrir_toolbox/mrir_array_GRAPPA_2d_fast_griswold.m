function dat_full = mrir_array_GRAPPA_2d_fast_griswold(dat, acs, prot)
%MRIR_ARRAY_GRAPPA_2D_FAST_GRISWOLD  wrapper around "fast_grappa3D"
%
% dat_full = mrir_array_GRAPPA_2d_fast_griswold(dat, acs, prot)


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/08
% $Id: mrir_array_GRAPPA_2d_fast_griswold.m,v 1.2 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  sig_fast = permute(squeeze(dat), [3, 2, 1, 4]);
  acs_fast = permute(squeeze(acs), [3, 2, 1, 4]);

  [dat_fast, ws_fast] = fast_grappa3D_2x2x3_sig(sig, acs, prot.lAccelFactPE, prot.lAccelFact3D);

  dat_full = permute(dat_fast, [3, 2, 1, 4]);

  dims = size(dat_full);

  % add back singleton dimensions to conform to ICE convention for ordering
  dat_full = reshape(raw, [dims(1), dims(2), dims(3), 1, 1, 1, 1, 1, dims(4)]);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_fast_griswold.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
