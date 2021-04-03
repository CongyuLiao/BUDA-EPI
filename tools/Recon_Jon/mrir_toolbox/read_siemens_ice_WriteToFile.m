function ima = read_siemens_ice_WriteToFile(imafile, filetype, varargin)
%READ_SIEMENS_ICE_WRITETOFILE
%
%  ima = read_siemens_ice_WriteToFile(imafile, filetype)
%
%  example:
%
%     ima_0001 = read_siemens_ice_WriteToFile('WriteToFile_0001.ima', 'short');
%

% NOTE: rumor has it that in the new ICE version, VB15, the offline
% simulator will generate DICOMs as if they were generated on the host. this
% long-overdue feature will make life much easier---and will render this
% function obsolete!   :)

% NOTE: this function was previously known as "read_siemens_ice_image.m".

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 10/05/2006
% $Id: read_siemens_ice_WriteToFile.m,v 1.1 2007/02/20 20:07:18 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  %------------------------------------------------------------------------%

  BASELINE = 'N4_VB13A_LATEST_20060607';

  if ( ~exist(imafile, 'file') ),
    errstr = sprintf('image file {%s} not found', imafile);
    error('!!! [%s]: %s', mfilename, errstr);
  end;

  fid = fopen(imafile, 'r', 'l');
  if ( fid == -1 ),
    errstr = sprintf('could not open image file {%s}', imafile);
    error('!!! [%s]: %s', mfilename, errstr);
  end;

  
  switch lower(filetype),
   case {'float32', 'float'},
    ima = fread(fid, 'float32');
   case {'short', 'uint16'},
    ima = fread(fid, 'uint16');  
   otherwise,
    errstr = sprintf('file type [%s] unrecognized', filetype);
    error('!!! [%s]: %s', mfilename, errstr);
  end;

  fclose(fid);

  if ( nargin == 2 ),
    [nRows, nCols] = deal(sqrt(length(ima)));
    if ( mod(nRows, 1) ~= 0 ),
      error('data is not a square matrix---specify "nRows" and "nCols"');
    end;
  end;

  if ( nargin == 3 ),
    dims = varargin{1};
    nRows = dims(1);
    nCols = dims(2);
  end;

  if ( nargin >= 4 ),
    nRows = varargin{1};
    nCols = varargin{2};
  end;

  if ( nRows*nCols ~= length(ima) ),
    warning('specified "nRows" and "nCols" does not match data in dump file');
    return;
  end;

  % dump runs first along first dimension (usually COL), then second
  % (usually LIN), so we transpose to get into [COL,LIN] into [rows,cols].
  ima = reshape(ima, nCols, nRows).';

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IDEA_ICE_DEV/ICE_TOOLS/MATLAB/read_siemens_ice_WriteToFile.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
