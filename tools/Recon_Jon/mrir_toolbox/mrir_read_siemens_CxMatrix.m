function cxmatrix = mrir_read_siemens_CxMatrix(filename, varargin)
% MRIR_READ_SIEMENS_CXMATRIX
%
% cxmatrix = mrir_read_siemens_CxMatrix(filename, nRows, nCols)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/20
% $Id: mrir_read_siemens_CxMatrix.m,v 1.1 2007/08/28 18:52:05 jonp Exp $
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
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_read_siemens_CxMatrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

