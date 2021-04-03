function img_uncomb = mrir_conventional_2d(raw, varargin);
%MRIR_CONVENTIONAL_2D  reconstructs conventional (cartesian) acquisitions
%
% img_uncomb = mrir_conventional_2d(raw);
%
% img_uncomb = mrir_conventional_2d(raw, prot);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id: mrir_conventional_2d.m,v 1.1 2007/09/13 19:46:19 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  prot = [];
  if ( nargin >= 2 ), prot = varargin{1}; end;
  if ( isempty(prot) ), prot = read_meas_prot__struct; end;

  if ( mrir_ice_dimensions(raw, 'par') > 1 ),
    warning(sprintf('input "%s" contains data along dimension "par"', inputname(1)));
  end;


  %==--------------------------------------------------------------------==%

  hyb_roft = mrir_iDFT_freqencode(raw);
  img_peft = mrir_iDFT_phasencode(hyb_roft, 'lin', prot.lPhaseEncodingLines);

  img_uncomb = mrir_image_crop(img_peft, prot.flReadoutOSFactor);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_conventional_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
