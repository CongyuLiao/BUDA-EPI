function [meas, FLAG__fixed] = mrir_array_channel_disable(meas, disabled_channels)
%MRIR_ARRAY_CHANNEL_DISABLE
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/15
% $Id: mrir_array_channel_disable.m,v 1.4 2011/03/28 04:14:45 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( isempty(disabled_channels) ),
    return;
  end;

  FLAG__fixed = 0;

  % assume that all meas.dat files have a "data" field
  if ( ~isfield(meas, 'data') ),
    error('no "data" field found');
  end;

  Ncha = mrir_ice_dimensions(meas.data, 'cha');

  cha = setdiff(1:Ncha, disabled_channels);



  for fieldnames = {'data', 'data_phascor1d', 'data_phascor2d', 'patrefscan', 'patrefscan_phascor', 'noiseadjscan', 'phasestabscan', 'phasestabtime', 'refphasestabscan', 'refphasestabtime', 'patrefphasestabscan', 'patrefphasestabtime'},
    fieldstr = char(fieldnames);
    if ( isfield(meas, fieldstr) ),

      measdata = getfield(meas, fieldstr);

      if ( mrir_ice_dimensions(measdata, 'cha') ~= Ncha ),
        disp(sprintf('==> [%s]: skipping field "%s" -- \t %d ~= %d', mfilename, ...
                     fieldstr, mrir_ice_dimensions(measdata, 'cha'), Ncha));
        continue;
      end;


      %                   1 2  3  4 5 6 7 8 9 0 1 2 3 4 5 6
      measdata = measdata(:,:,cha,:,:,:,:,:,:,:,:,:,:,:,:,:);

      meas = setfield(meas, fieldstr, measdata);

      FLAG__fixed = 1;

    end;  %% if
  end;  %% for


  if ( FLAG__fixed ),
    disp(sprintf('==> [%s]: removed the following disabled channels from data: %s', mfilename, num2str(disabled_channels, '%02d ')));
    meas.evp.NChaMeas = length(cha);
  end;



  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_channel_disable.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
