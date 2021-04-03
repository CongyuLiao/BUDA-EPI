function D = mrir_iDFT2_matrix(N1, N2)
%MRIR_IDFT2_MATRIX  efficient calculation of *unitary-scaled* 2D DFT matrix
%
% D = mrir_iDFT2_matrix(N1, N2)
%
%  
%  usage: 
%    
%     x = ifft2(k) * (sqrt(N1)*sqrt(N2));
%
%  is equivalent to
%
%     x = reshape(D*k(:), [N1, N2])

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/jan/11
% $Id: mrir_iDFT2_matrix.m,v 1.2 2009/01/30 06:26:48 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  
  % 'unitary' version of matrix, so adjoint(D) = hermitian(D)
  D1 = conj(dftmtx(N1))/sqrt(N1); D2 = conj(dftmtx(N2))/sqrt(N2);
  D = kron(D2, D1);
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT2_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
