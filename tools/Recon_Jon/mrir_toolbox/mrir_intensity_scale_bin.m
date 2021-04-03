function img_scale_bin = mrir_intensity_scale_bin(img, scale)
%MRIR_INTENSITY_SCALE_BIN
%
% img_scale_bin = mrir_intensity_scale_bin(img, scale)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/23
% $Id: mrir_intensity_scale_bin.m,v 1.1 2007/01/23 23:06:30 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( ~exist('scale', 'var') ),
    % 4095
    scale = (2^12) - 1;
  end;
  
  img_scale = ( img - min(img(:)) ) ./ ( max(img(:)) - min(img(:)) ) * scale;
  img_scale_bin = uint16(round( img_scale ));

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_intensity_scale_bin.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
