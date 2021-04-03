function [meas, FLAG__fixed] = mrir_hack__fix_7T_channels(meas, varargin)
%MRIR_HACK__FIX_7T_CHANNELS  remove empty channel 1 in 32-channel coil data
%
% meas = mrir_hack__fix_7T_channels(meas)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/jan/07
% $Id: mrir_hack__fix_7T_channels.m,v 1.4 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__fixed = 0;

  if ( ~isfield(meas.prot, 'ManufacturersModelName') ),
    % quietly abort
    return;
  end;

  if ( isempty(meas.evp.RawCha) ),
    c = sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(meas.data, 1), 2), 4), 5), 6), 7), 8), 9), 10), 11), 12), 13), 14), 15), 16);
    meas.evp.RawCha = length(find(c));
  end;
  
  
  if ( strcmp(meas.prot.ManufacturersModelName, 'Investigational_Device_7T') ),
    if ( (meas.evp.RawCha < 32) && (mrir_ice_dimensions(meas.data, 'cha') == 32) ),

      % allow user to specify additional channels to disable
      if ( nargin >= 2 ),
        disabled_channels = varargin{1};
      else,
	disabled_channels = [];
      end;

      % disable first channel and any that user specified
      [meas, FLAG__fixed] = mrir_array_channel_disable(meas, [1; disabled_channels(:)]);


    end;  %% if
  else,
    warning(sprintf('meas data not from 7T: "%s"', meas.prot.ManufacturersModelName));
  end;  %% if


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_hack__fix_7T_channels.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
