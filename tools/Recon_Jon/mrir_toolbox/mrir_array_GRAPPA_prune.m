function [datprune, acsprune] = mrir_array_GRAPPA_prune(dat, acs, evp)
%MRIR_ARRAY_GRAPPA_PRUNE  extract data and ACS lines from sparse arrays
%
% [datprune, acsprune] = mrir_array_GRAPPA_prune(dat, acs, evp)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/06
% $Id: mrir_array_GRAPPA_prune.m,v 1.2 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

VERSION = '$Revision: 1.2 $';
if ( nargin == 0 ), help(mfilename); return; end;


%==--------------------------------------------------------------------==%

[datlines1, acslines1] = mrir_array_GRAPPA_indices(evp.NAFLin, ...
    evp.NFirstLin, evp.NLinMeas, ...
    evp.NFirstRefLin, evp.NRefLin);

[datlines2, acslines2] = mrir_array_GRAPPA_indices(evp.NAFPar, ...
    evp.NFirstPar, evp.NParMeas, ...
    evp.NFirstRefPar, evp.NRefPar);

if ( ~isempty(dat) ),
    %              1         2 3 4 5 6 7 8         9 0 1 2 3 4 5 6
    datprune = dat(:,datlines1,:,:,:,:,:,:,datlines2,:,:,:,:,:,:,:);
    mrir_array_GRAPPA_check_density(sum(datprune, 8));
    datCheck = 1;
else,
    datprune = [];
    datCheck = 0;
end;


if ( ~isempty(acs) ),
    %              1         2 3 4 5 6 7 8         9 0 1 2 3 4 5 6
    acsprune = acs(:,acslines1,:,:,:,:,:,:,acslines2,:,:,:,:,:,:,:);
    %acsprune = single(acs(:,acslines1,:,:,:,:,:,:,acslines2,:,:,:,:,:,:,:));
    mrir_array_GRAPPA_check_density(sum(acsprune, 8));
else,
    acsprune = [];
end;


%==--------------------------------------------------------------------==%
if datCheck == 1
    if ( mrir_ice_dimensions(datprune, 'lin') ~= evp.RawLin )
        error('number of data lines = %d, expected data lines = [%d]', ...
            mrir_ice_dimensions(datprune, 'lin'), evp.RawLin);
    end;
    
    % another nice inconsistency
    
    if ( mrir_ice_dimensions(datprune, 'par') ~= evp.NParMeas ),
        warning('number of data partitions = %d, expected data partitions = [%d]', ...
            mrir_ice_dimensions(datprune, 'par'), evp.NParMeas);
    end;
    
end
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_prune.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
