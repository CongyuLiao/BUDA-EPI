function phases_receivers_report = mrir_array_modemix_calibrate_import(dirname)
%MRIR_ARRAY_MODEMIX_CALIBRATE_IMPORT  
%
% phases_receivers_report = mrir_array_modemix_calibrate_import(dirname)
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/oct/26
% $Id: mrir_array_modemix_calibrate_import.m,v 1.1 2008/10/27 00:54:36 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  keyword = 'ID_RFLinMeanPhase';
  
      
  cmd = sprintf('grep -H -A 1 %s %s/*.xml | grep DoubleVariable', keyword, dirname);
  [status, result] = system(cmd);
  
  v = regexp(result, '"(?<value>-?\d+\.*\d*)"', 'names');
  
  if ( isempty(v) ),
    error('==> [%s]: modulator/receiver phases not found in specified directory', mfilename);
  end;
  
  for ind = 1:length(v),
    phases_receivers_report(ind) = str2num(v(ind).value);
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_modemix_calibrate_import.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
