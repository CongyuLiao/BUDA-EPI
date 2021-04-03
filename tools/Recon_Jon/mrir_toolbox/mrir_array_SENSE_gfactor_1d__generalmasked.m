function gfactor = mrir_array_SENSE_gfactor_1d__general(receive, covmtx, R)
%MRIR_ARRAY_SENSE_GFACTOR_1D__GENERAL
%
% gfactor = mrir_array_SENSE_gfactor_1d__general(receive, covmtx, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id: mrir_array_SENSE_gfactor_1d__general.m,v 1.2 2007/05/14 02:36:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(receive,  1); % frequency encoded
  Nlin = size(receive,  2); % phase encoded
  Ncha = size(receive,  3);
  Nslc = size(receive, 10); % slices

  % whitening operator for pre-whitening data
  W = mrir_array_whitening_operator(covmtx, 'svd');

  receive_whitened = reshape(receive, [], Ncha) * W.';
  receive = reshape(receive_whitened, size(receive));

  % image folding operator, F: ceil(Nlin/R) x Nlin
  F = mrir_array_accelerated_folding_operator(Nlin, R);
  
  % preallocate
  gfactor = zeros(Ncol, Nlin);

  for slc = 1:Nslc,
  
  for ii = 1:Ncol,
     
      % receive matrix, S: Nlin x Ncha
      S = squeeze(receive(ii, :, :, 1,1,1,1,1,1, slc));

      index = sum(abs(S),2) > 0; %mask out area where there is no sensitivity
      %errcovmtx_inv = zeros(Nlin, Nlin);
      if sum(index) ~= 0
          errcovmtx_inv = zeros(sum(index), sum(index));
          for jj = 1:ceil(Nlin/R),

              % strip of row of folding operator corresponding to this line
              %a = F(jj,:);
              a = F(jj,index);

              % reception encoding matrix (E: Ncha x Nlin) computes the weighted sum
              % of all "true" pixels in this line for the observed pixel value at
              % each coil. each line is weighted first by the receive profile of a
              % coil, then the weighted pixels are summed according to the aliasing
              % for this row. "a" is sparse, and the number of nonzero elements is
              % greater than or equal to R.
              E = (diag(a) * S(index,:)).';

              % since in general the pixels contributing to an observation are not
              % unique to the observation, the errors must be summed.
              errcovmtx_inv = errcovmtx_inv + (E'*E);

          end;


          %    if ( cond(errcovmtx_inv) > 1/sqrt(eps) ),
          %      error('the error covariance matrix is singular!');
          %    end;
          errcovmtx = inv(errcovmtx_inv);

          % store g-factor for all lines in this image column
          gfactor(ii, index, 1, 1,1,1,1,1,1, slc) = sqrt(abs( diag(errcovmtx) .* diag(errcovmtx_inv) ));
      end
  end;

  end;
  
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__general.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
