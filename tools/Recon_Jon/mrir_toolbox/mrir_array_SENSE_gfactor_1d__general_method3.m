function gfactor = mrir_array_SENSE_gfactor_1d__general_method3(receive, covmtx, R)
% combine "method 2" with "special" case one 

%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(receive,  1); % frequency encoded
  Nlin = size(receive,  2); % phase encoded
  Ncha = size(receive,  3);
  Nslc = size(receive, 10); % slices

  gfactor = zeros(Ncol, Nlin, 1, 1,1,1,1,1,1, Nslc);
  
  % image folding operator, F: ceil(Nlin/R) x Nlin
  F = mrir_array_accelerated_folding_operator(Nlin, R);
 if  (sum (sum( (abs(F) > eps*10) ,1)  > 1) > 0) %i.e. an unaliased pixel can be from more than 1 aliased pixel
      % preallocate
      covmtx_h = kron(eye(length(1:R:Nlin)),covmtx);
      covmtxinv_h = inv(covmtx_h);

      for slc = 1:Nslc,

          for ii = 1:Ncol,

              % receive matrix, S: Nlin x Ncha
              S = squeeze(receive(ii, :, :, 1,1,1,1,1,1, slc));
              index = sum(abs(S),2) > 0; %mask out area where there is no sensitivity
              if sum(index) ~= 0

                  E = [];
                  for jj = 1:size(F,1) % floor(Nlin/R),

                      % strip of row of folding operator corresponding to this line
                      a = F(jj,index);

                      % reception encoding matrix (E: Ncha x Nlin) computes the weighted sum
                      % of all "true" pixels in this line for the observed pixel value at
                      % each coil. each line is weighted first by the receive profile of a
                      % coil, then the weighted pixels are summed according to the aliasing
                      % for this row. "a" is sparse, and the number of nonzero elements is
                      % greater than or equal to R.
                      E_current = (diag(a) * S(index,:)).';

                      E  = [ E; E_current];
                  end;

                  errcovmtx_inv = E' * covmtxinv_h * E;
                  errcovmtx = inv(errcovmtx_inv);
                  
    
                  % store g-factor for all lines in this image column
                  gfactor(ii, index, 1, 1,1,1,1,1,1, slc) = sqrt(abs( diag(errcovmtx) .* diag(errcovmtx_inv) ));
                  %I(ii, index, 1, 1,1,1,1,1,1, slc) = (errcovmtx) * E' * covmtxinv_h * reshape(I_Collapsed(ii,index,:, 1,1,1,1,1,1, slc),[],1);
              end
          end;
      end;
  else
      covmtxinv = inv(covmtx);

      aliasmap = mrir_array_accelerated_aliasmap(Nlin, R);

      for slc = 1:Nslc,

          sens = receive(:,:,:, 1,1,1,1,1,1, slc);

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
                      %I(ii, aliased_pixels_ind(index), 1, 1,1,1,1,1,1, slc) = (errcovmtx) * E' * covmtxinv * reshape(I_Collapsed(ii,aliased_pixels_ind(index),:, 1,1,1,1,1,1, slc),[],1);
                  end
              end;
          end;
      end;
  end

     
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__general.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
