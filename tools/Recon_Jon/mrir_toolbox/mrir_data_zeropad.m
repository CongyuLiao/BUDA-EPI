function padded = mrir_data_zeropad(meas, dim, finalsize)
%MRIR_DATA_ZEROPAD  zero pad arrays to fill in skipped lines/partitions
%
% padded = mrir_data_zeropad(meas, dim, finalsize)
%
%
% example:  
%
%   meas.patrefscan = mrir_data_zeropad(meas.patrefscan, 2, meas.prot.iNoOfFourierLines);

% TODO: generalize for cases where input is a struct with multiple fields
% needing padding and where input struct also contains the prot struct.
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/aug/28
% $Id: mrir_data_zeropad.m,v 1.1 2007/08/28 17:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  padamount = repmat(0, 1, ndims(meas));
  padamount(dim) = finalsize - size(meas, dim);

  padded = padarray(meas, padamount, 'post');


  return;
  
  
  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_data_zeropad.m,v $
  
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
