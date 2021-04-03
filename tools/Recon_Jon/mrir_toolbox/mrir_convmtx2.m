function C = mrir_convmtx2(f, N1, N2)
%MRIR_CONVMTX2  wrapper around "convmtx2" for toeplitz matrices (conv2 'same')
%
% C = mrir_convmtx2(f, N1, N2)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/feb/16
% $Id: mrir_convmtx2.m,v 1.1 2009/02/16 06:06:29 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  m1 = size(f, 1);
  m2 = size(f, 2);

  F = convmtx2(f, N1, N2);


  L1 = N1+m1-1;
  L2 = N2+m2-1;

  delta1 = (m1-1)/2;
  delta2 = (m2-1)/2;


  [row, col] = meshgrid( 1+delta2 : L2-delta2, 1+delta1 : L1-delta1 );

  ind = sub2ind([L1 L2], col(:), row(:));

  C = F(ind,:);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_convmtx2.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
