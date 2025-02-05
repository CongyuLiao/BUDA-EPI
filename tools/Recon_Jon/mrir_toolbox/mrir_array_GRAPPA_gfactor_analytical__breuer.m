function [g, g_ind, n, n_ind, covmtx] = mrir_array_GRAPPA_gfactor_analytical__breuer(combine_weights, G, i_noisecov, NColMeas, NLinMeas, NParMeas)
%MRIR_ARRAY_GRAPPA_GFACTOR_ANALYTICAL__BREUER
%
% [g, g_ind] = mrir_array_GRAPPA_gfactor_analytical__breuer(combine_weights, G, i_noisecov, NColMeas, NLinMeas, NParMeas)

% Breuer et al. (2008), "A general formulation for quantitative g-factor
% calculation in GRAPPA reconstructions [Abstract]"; Proc. ISMRM 16:10.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/13
% $Id: mrir_array_GRAPPA_gfactor_analytical__breuer.m,v 1.3 2009/03/03 02:24:00 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  NCha = unique(size(i_noisecov));

  if ( ~isstruct(G) ),

    disp('using pre-computed kernel');
    K = G;
  else,
    k = mrir_array_GRAPPA_conv_kernel(G);
    K = mrir_image_crop(mrir_iDFT(mrir_iDFT(k, 1, NColMeas), 2, NLinMeas), 1);
  end;

  flReadoutOSFactor = mrir_ice_dimensions(K, 'col') ./ mrir_ice_dimensions(combine_weights, 'col');
  K = mrir_image_crop(K, flReadoutOSFactor);


  i_noisestd = sqrt(diag(i_noisecov));

  for par = 1:NParMeas,
    for lin = 1:NLinMeas,
      for col = 1:NColMeas/flReadoutOSFactor,

        % matrix W: [Ctrg x Csrc]
        W = squeeze(K(col, lin, par, :, :)).';

        covmtx(col, lin, :,:) = conj(W) * i_noisecov * W.' ;
        n_ind(col, lin, :) = diag(sqrt(abs( conj(W) * i_noisecov * W.' )));
        g_ind(col, lin, :) = squeeze(n_ind(col, lin, :)) ./ i_noisestd;

        % vector p: [1 x Ctrg]        1    2  3  4 5 6 7 8    9
        p = squeeze(combine_weights(col, lin, :, 1,1,1,1,1, par))';

        n(col, lin) = sqrt(abs( p * W * i_noisecov * W' * p' ));
        g(col, lin) = n(col, lin) ./ sqrt( abs(p * i_noisecov * p') );

      end;
    end;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_gfactor_analytical__breuer.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
