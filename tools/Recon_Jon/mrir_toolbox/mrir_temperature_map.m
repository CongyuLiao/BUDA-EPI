function img_temp = mrir_temperature_map(img01, img02, B0, TE)
%MRIR_TEMPERATURE_MAP  compute temperature change between two image acquisitions from phase
%
% img_temp = mrir_temperature_map(img01, img02, B0, TE)

% references:
%
%   Shapiro EM, Borthakur A, Shapiro MJ, Reddy R, Leigh JS.
%   Fast MRI of RF heating via phase difference mapping.
%   Magn Reson Med. 2002 Mar;47(3):492-8.
%   PMID: 11870836
%   
%   Shapiro EM, Borthakur A, Reddy R.
%   MR imaging of RF heating using a paramagnetic doped agarose phantom.
%   MAGMA. 2000 Jun;10(2):114-21.
%   PMID: 10873201
%   
%   MacFall JR, Prescott DM, Charles HC, Samulski TV.
%   1H MRI phase thermometry in vivo in canine brain, muscle, and tumor tissue.
%   Med Phys. 1996 Oct;23(10):1775-82.
%   PMID: 8946373
%   
%   Ishihara Y, Calderon A, Watanabe H, Okamoto K, Suzuki Y, Kuroda K, Suzuki Y.
%   A precise and fast temperature mapping using water proton chemical shift.
%   Magn Reson Med. 1995 Dec;34(6):814-23.
%   PMID: 8598808
%  
%   Chen CN, Hoult DI.
%   The visualization of RF probe electric fields.
%   Magn Reson Med. 1993 Mar;29(3):386-90.
%   PMID: 8450747
    
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/11
% $Id: mrir_temperature_map.m,v 1.3 2008/08/25 05:30:31 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  %% TODO (1): allow passing parameters as input
  
  % K: "temperature dependent chemical shift coefficient in ppm/degC"
  % (for water, K is approximately 0.01 ppm, [shapiro2002fast])
  K  = -0.011e-6;

  gamma = 42.576e6;  % Hz/T

%  B0 = 7.0;
%  TE = 8e-3;


  %% TODO (2): determine absolute phase (offset from B0) to get absolute temperature


  %==--------------------------------------------------------------------==%
  % phase unwrapping through calls to FSL's "prelude"

  %% TODO (3): read in mask used by "prelude" (crashes for some reason)

  tmpfile = fullfile(tempdir, [mfilename, '__', num2str(now)], 'CREATED');
  tmpdir = fileparts(tmpfile);
  [status, result] = system(sprintf('mkdir -p %s', tmpdir));
  if ( status ), error(result); end;
  [status, result] = system(sprintf('touch %s', tmpfile));
  if ( status ), error(result); end;


  save_nii(make_nii(squeeze(abs(img01))),   sprintf('%s/abs01.nii', tmpdir));
  save_nii(make_nii(squeeze(angle(img01))), sprintf('%s/phz01.nii', tmpdir));

  save_nii(make_nii(squeeze(abs(img02))),   sprintf('%s/abs02.nii', tmpdir));
  save_nii(make_nii(squeeze(angle(img02))), sprintf('%s/phz02.nii', tmpdir));


  prelude = '/usr/pubsw/packages/fsl/current/bin/prelude';

  [status, result] = system(sprintf('%s -a %s/abs01.nii -p %s/phz01.nii -s -u %s/phz01_unwrap.nii', prelude, tmpdir, tmpdir, tmpdir));
  if ( status ), error(result); end;
  [status, result] = system(sprintf('gunzip -c %s/phz01_unwrap.nii.gz > %s/phz01_unwrap.nii', tmpdir, tmpdir));
  if ( status ), error(result); end;

  [status, result] = system(sprintf('%s -a %s/abs02.nii -p %s/phz02.nii -s -u %s/phz02_unwrap.nii', prelude, tmpdir, tmpdir, tmpdir));
  if ( status ), error(result); end;
  [status, result] = system(sprintf('gunzip -c %s/phz02_unwrap.nii.gz > %s/phz02_unwrap.nii', tmpdir, tmpdir));
  if ( status ), error(result); end;

  phz01 = getfield(load_nii(sprintf('%s/phz01_unwrap.nii', tmpdir)), 'img');
  phz02 = getfield(load_nii(sprintf('%s/phz02_unwrap.nii', tmpdir)), 'img');


  %==--------------------------------------------------------------------==%

  mask_bin = ( (phz01.*phz02) ~= 0 );


  for slc = 1:size(mask_bin, 3),
    mask_bin(:,:,slc) = bwmorph(mask_bin(:,:,slc), 'dilate', 1);
  end;

  mask_ind = find( mask_bin );


  dphz = phz02 - phz01;

  %% TODO (4): remove this hack
  
  % this is a hack to correct strange problem where some slices have 2PI jump
  for slc = 1:size(dphz, 3),
    dphz_slc = dphz(:,:,slc);
    if ( mean( dphz_slc(find(dphz_slc)) ) < 0 )
      dphz_slc(find(dphz_slc)) = dphz_slc(find(dphz_slc)) + 2 * pi;
%      disp(sprintf('found %02d', slc));
    end;
    dphzC(:,:,slc) = dphz_slc;
  end;

  % phase-temperature equation [see macfall1996h1]
  img_temp = rad2deg(dphzC) / 360 / K / B0 / gamma / TE;

  
  %% TODO (5): remove this hack
  
  % hack to get temperature into reasonable range
  offset = median( img_temp(find(img_temp)) );
  img_temp(find(img_temp)) = img_temp(find(img_temp)) - offset;


  return;

  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_temperature_map.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
