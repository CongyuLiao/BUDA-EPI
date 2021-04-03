function birdcage_modes = mrir_array_birdcage(strips)
%MRIR_ARRAY_BIRDCAGE
%
% birdcage_modes = mrir_array_birdcage(strips)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/18
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = mrir_ice_dimensions(strips, 'col');
  Nlin = mrir_ice_dimensions(strips, 'lin');
  Ncha = mrir_ice_dimensions(strips, 'cha');

  delta_phase = 360 / Ncha;


  mode01_phases = -deg2rad([ delta_phase : delta_phase : 360 ]);

  birdcage_modes = zeros(Ncol, Nlin, Ncha);

  for mode = 1:Ncha,

    birdcage_chan = zeros(Ncol, Nlin);
    for cha = 1:Ncha,
      birdcage_chan(:,:,cha) = strips(:,:,cha) * exp( i*mode01_phases(cha) * mode );
    end;

    birdcage_modes(:,:,mode) = mean(birdcage_chan, 3);

  end;


  return;

  

  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
