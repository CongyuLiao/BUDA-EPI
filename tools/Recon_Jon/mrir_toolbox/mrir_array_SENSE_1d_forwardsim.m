function img_redu = mrir_array_SENSE_1d_forwardsim(img_full, sensitivity, R)
%MRIR_ARRAY_SENSE_1D_FORWARDSIM  simulate sensitivity encoding (for debugging)
%
% img_redu = mrir_array_SENSE_1d_forwardsim(img_full, sensitivity, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jul/07
% $Id: mrir_array_SENSE_1d_forwardsim.m,v 1.1 2008/05/31 02:59:24 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  img_redu = zeros(Ncol, Nlin/R);

  for jj = 1:Nlin/R,

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);

    for ii = 1:Ncol,
      for ll = 1:Ncha,
        img_redu(ii, jj, ll) = sensitivity(ii, aliasmap(jj,:), ll) * ...
            img_full(ii, aliasmap(jj,:)).';
      end;
    end;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_1d_forwardsim.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
