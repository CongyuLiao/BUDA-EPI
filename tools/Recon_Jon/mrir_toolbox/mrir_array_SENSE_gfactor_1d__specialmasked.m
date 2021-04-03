function gfactor = mrir_array_SENSE_gfactor_1d__specialmasked(sensitivity, covmtx, R)
%MRIR_ARRAY_GFACTOR_1D__SPECIAL
%
% gfactor = mrir_array_Gfactor_1d__special(sensitivity, covmtx, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id: mrir_array_SENSE_gfactor_1d__special.m,v 1.1 2007/05/05 03:35:28 jonp Exp $

% Kawin Setsompop : do gfactor with account of masking


%**************************************************************************%

% I = mrir_array_combine(sensitivity,0);
% sensitivity( repmat( I < mean(I(:))/5,[1 1 size(sensitivity,3)] ) ) = 0;

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity,  1); % frequency encoded
  Nlin = size(sensitivity,  2); % phase encoded
  Ncha = size(sensitivity,  3);
  Nslc = size(sensitivity, 10);
  
  covmtxinv = inv(covmtx);

  aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

  % preallocate
  gfactor = zeros(Ncol, Nlin, 1, 1,1,1,1,1,1, Nslc);

  for slc = 1:Nslc,

    sens = sensitivity(:,:,:, 1,1,1,1,1,1, slc);
    
  for jj = 1:Nlin/R,

    % pull out the indices of the aliased pixels for this position from lookup table
    aliased_pixels_ind = aliasmap(jj, 1:R);
  
    for ii = 1:Ncol,
        
        index = squeeze(sum( abs(  sens(ii,aliased_pixels_ind,:) ) ,3) ) > 0; %only use pixel with non zero sensitivity
        if sum(index) ~= 0
            % encoding matrix, E: Ncha x sum(index) where sum(index) is the number of overlapping pixels
            E = permute(reshape( sens(ii, aliased_pixels_ind(index), :), sum(index), Ncha ), [2 1]);

            errcovmtx_inv = E' * covmtxinv * E;

            %      if ( cond(errcovmtx_inv) > 1/sqrt(eps) ),
            %        error('the error covariance matrix is singular!');
            %      end;
            errcovmtx = inv(errcovmtx_inv);

            gfactor(ii, aliased_pixels_ind(index), 1, 1,1,1,1,1, slc) = sqrt(abs( diag(errcovmtx) .* diag(errcovmtx_inv) ));
        end
    end;
  end;
  end;
  
  return;


  
  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__special.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: