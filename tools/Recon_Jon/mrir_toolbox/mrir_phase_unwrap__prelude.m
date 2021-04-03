function phz_unwrap = mrir_phase_unwrap__prelude(img)
%MRIR_PHASE_UNWRAP__PRELUDE
%
% phz_unwrap = mrir_phase_unwrap__prelude(phz)

% this code began as a part of the function "mrir_temperature_map.m".

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/21
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  % TODO: use global variable PRELUDE as a matlab environment-like variable
  prelude = '/usr/pubsw/packages/fsl/current/bin/prelude';


  %==--------------------------------------------------------------------==%
  % phase unwrapping through calls to FSL's "prelude"

  tmpdir = mrir_sysutil__tempdir(mfilename);


  save_nii(make_nii(squeeze(  abs(img))), sprintf('%s/abs.nii', tmpdir));
  save_nii(make_nii(squeeze(angle(img))), sprintf('%s/phz.nii', tmpdir));


  [status, result] = system(sprintf(...
      '%s -v -a %s/abs.nii -p %s/phz.nii -s -u %s/phz_unwrap.nii', ...
       prelude, tmpdir,       tmpdir,          tmpdir));

  if ( status ), error(result); end;


  [status, result] = system(sprintf(...
      'gunzip -c %s/phz_unwrap.nii.gz > %s/phz_unwrap.nii', ...
                 tmpdir,                tmpdir));

  if ( status ), error(result); end;

  phz_unwrap = getfield(load_nii(sprintf('%s/phz_unwrap.nii', tmpdir)), 'img');


  %% TODO: read in mask used by "prelude" (crashes for some reason)


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
