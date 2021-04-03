function header = read_meas_dat__header(filename)
%READ_MEAS_DAT__HEADER  return header from "meas.dat" file as ASCII string
%
% header = read_meas_dat__header(filename)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/27
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [fp, errstr] = fopen(filename, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;

  % determine size (in bytes) of ascii header files stored in the 'meas.dat'
  % format (i.e., "Config_.evp", "Dicom_.evp", etc) to skip over them all.
  % [note: this works for VB11A also---the first integer is 32, which is
  % the number of bytes to skip at the beginning!]
  data_start = fread(fp, 1, 'uint32');

  % read header into one string for parsing
  header = fscanf(fp, '%c', data_start-4);


  fclose(fp);


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
