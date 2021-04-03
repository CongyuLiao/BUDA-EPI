function mrir_array_modemix_attenuators_write(filename, header, modemixmatrix_atten_steps_dB)
%MRIR_ARRAY_MODEMIX_ATTENUATORS_WRITE
%
% mrir_array_modemix_attenuators_write(filename, header, modemixmatrix_atten_steps_dB)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/oct/26
% $Id: mrir_array_modemix_attenuators_write.m,v 1.1 2008/10/26 22:32:32 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  fp = fopen(filename, 'w');

  fprintf(fp, '%s  [%s]\n', header, datestr(now));

  for board = 1:16,
    fprintf(fp, '%02d    ', board);
    fprintf(fp, '%5.1f  ', modemixmatrix_atten_steps_dB(:,board));
    fprintf(fp, '\n');
  end;

  fclose(fp);


  return;

  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_modemix_attenuators_write.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
