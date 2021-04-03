function g = mrir_array_SENSE_gfactor_1d__arne(img, sens, covmtx, R, reg1, reg2)
%MRIR_ARRAY_SENSE_GFACTOR_1D__ARNE
%
% g = mrir_array_SENSE_gfactor_1d__arne(img, sens, covmtx, R, reg1, reg2)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/18
% $Id: mrir_array_SENSE_gfactor_1d__arne.m,v 1.1 2008/12/18 21:30:00 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  NLin = mrir_ice_dimensions(img, 'lin');


  for lin = 1:NLin,
    SNR_LOSS(lin, :) = arne_gfactor_SNR_code__jrp(squeeze(img(lin,:,:)), squeeze(sens(lin,:,:)), R, covmtx, reg1, reg2);
  end;

  g = 1./SNR_LOSS;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__arne.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
