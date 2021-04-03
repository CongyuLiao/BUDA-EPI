function img_uncomb = mrir_conventional_3d(raw, varargin)
%MRIR_CONVENTIONAL_3D  reconstructs conventional (cartesian) acquisitions
%
% img_uncomb = mrir_conventional_3d(raw);
%
% img_uncomb = mrir_conventional_3d(raw, prot);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id: mrir_conventional_3d.m,v 1.1 2007/09/13 19:46:47 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  prot = [];
  if ( nargin >= 2 ), prot = varargin{1}; end;
  if ( isempty(prot) ), prot = read_meas_prot__struct; end;

  if ( mrir_ice_dimensions(raw, 'slc') > 1 ),
    warning(sprintf('input "%s" contains data along dimension "slc"', inputname(1)));
  end;

  varstruct = whos('raw');
  if ( (varstruct.bytes / 2^30) > 8 ),
    warning('input exceeds 8 GB, FFT overflow may occur!');
  end;
  

  %==--------------------------------------------------------------------==%

  hyb_roft = mrir_iDFT_freqencode(raw);
  img_peft = mrir_iDFT_phasencode(hyb_roft, 'lin', prot.lPhaseEncodingLines);
  img_paft = mrir_iDFT_phasencode(img_peft, 'par', prot.lPartitions);

  img_uncomb = mrir_image_crop(img_paft, prot.flReadoutOSFactor);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_conventional_3d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
