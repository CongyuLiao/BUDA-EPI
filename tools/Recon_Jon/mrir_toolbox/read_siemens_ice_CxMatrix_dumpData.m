function cxmatrix = read_siemens_ice_CxMatrix_dumpData(filename, varargin)
% READ_SIEMENS_ICE_CXMATRIX_DUMPDATA
%
% cxmatrix = read_siemens_ice_CxMatrix_dumpData(filename, nRows, nCols)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/20
% $Id: read_siemens_ice_CxMatrix_dumpData.m,v 1.1 2007/07/17 23:50:13 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  fp = fopen(filename);
  dump = fread(fp, 'float');

  if ( nargin == 1 ),
    [nRows, nCols] = deal(sqrt(length(dump)/2));
    if ( mod(nRows, 1) ~= 0 ),
      error('data is not a square matrix---specify "nRows" and "nCols"');
    end;
  end;

  if ( nargin == 2 ),
    dims = varargin{1};
    nRows = dims(1);
    nCols = dims(2);
  end;

  if ( nargin >= 3 ),
    nRows = varargin{1};
    nCols = varargin{2};
  end;

  if ( nRows*nCols ~= length(dump)/2 ),
    error('specified "nRows" and "nCols" does not match data in dump file');
  end;

  % dump runs first along columns, then rows
  cxmatrix = complex(reshape(dump(1:2:end), nCols, nRows), reshape(dump(2:2:end), ...
                                                  nCols, nRows)).';

  fclose(fp);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_siemens_ice_CxMatrix_dumpData.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

