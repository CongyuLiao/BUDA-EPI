function dat_pad = mrir_image_pad_meas(dat, evp)
%MRIR_IMAGE_PAD_MEAS  pad arrays missing zeros at ends
%
% dat_pad = mrir_image_pad_meas(dat, evp)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/feb/14
% $Id: mrir_image_pad_meas.m,v 1.2 2009/02/15 00:21:25 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  NLin = mrir_ice_dimensions(dat, 'lin');
  NPar = mrir_ice_dimensions(dat, 'par');

  if ( NLin < evp.NLinMeas ),

    warning('zero-padding k-space data to desired number of lines.');
    dat_pad = padarray(dat, [0,(evp.NLinMeas-NLin),0,0,0,0,0,0], 0, 'post');

  else,
    dat_pad = dat;
  end;


  if ( NPar < evp.NParMeas ),

    warning('zero-padding k-space data to desired number of partitions.');
    dat_pad = padarray(dat_pad, [0,0,0,0,0,0,0,(evp.NParMeas-NPar)], 0, 'post');

  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_pad_meas.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
