function C = mrir_cconvmtx2(f, N1, N2)
%MRIR_CCONVMTX2  wrapper around "cconv2" for circulant matrices
%
% C = mrir_cconvmtx2(f, N1, N2)

% (based on "cconvmtx2.m" by W. Clem Karl)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/mar/10
% $Id: mrir_convmtx2.m,v 1.1 2009/02/16 06:06:29 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  m1 = size(f, 1);
  m2 = size(f, 2);

  C = sparse([], [], [], N1*N2, N1*N2, N1*N2*m1*m2);

  for jj = 1:N2,
    for ii = 1:N1,
      ej = zeros(N1, N2);
      ej(ii, jj) = 1;
      tmp = cconv2(f, ej, N1, N2);
      C(:, (jj-1)*N1+ii) = tmp(:);
    end;
  end;



  return;


%**************************************************************************%
function z = cconv2(x, y, N1, N2)


  z = ifft2( fft2(x, N1, N2) .* fft2(y, N1, N2) );

  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_convmtx2.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
