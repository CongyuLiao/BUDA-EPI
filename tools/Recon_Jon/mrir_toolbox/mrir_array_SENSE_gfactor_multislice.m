function gfactor = mrir_array_SENSE_gfactor_multislice(receive, covmtx, R)

%Kawin Setsompop 23 June 09
%
% calc gfactor for multislice case, only do one slice group at a time. 

%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(receive,  1); % frequency encoded
  Nlin = size(receive,  2); % phase encoded
  Ncha = size(receive,  3);
  Nslc = size(receive,  10); % slices

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
                  end
              end;
          end;
      end;
 end
  
% % % %  ??????
% % % %      covmtxinv = inv(covmtx);
% % % %     for jj = 1:size(F,1) %floor(Nlin/R)
% % % %         aliased_index = find(F(jj,:) > eps*10);
% % % %         for ii = 1:Ncol,
% % % %             % receive matrix, S: Nlin x Ncha
% % % %             S = receive(ii, aliased_index, :,:); % 1 x R_Inplane x Ncoils x Nslices
% % % %             index = sum(abs(S),3) > eps*10; %only involve pixels with non zero sensitivity
% % % %             if sum(index(:)) ~= 0 %i.e. if not all zero
% % % %                 E = [];
% % % %                 for SlcCount = 1:R_slc
% % % %                      E_current = permute( reshape( receive(ii, aliased_index(index(:,:,:,SlcCount)), :, SlcCount) , sum(index(:,:,:,SlcCount)), Ncha ), [2 1]);
% % % %                      E = [E, E_current]; %E: Ncha x #AliasedPixels
% % % %                 end
% % % %                 %noise-weighted pseudoinverse: Nalias x Nchan
% % % %                 U = inv( (E' * covmtxinv * E) ) * E' * covmtxinv;
% % % % 
% % % %                 % observed reduced-FOV image intensities: Nchan x 1
% % % %                 a = squeeze(img_redu(ii, jj, :));
% % % % 
% % % %                 % SENSE estimate: Nalias x 1
% % % %                 v = U * a;
% % % % 
% % % %                 StartIndex = 1;
% % % %                 EndIndex = sum(index(:,:,:,1));
% % % %                 for SlcCount = 1:R_slc
% % % %                     % stuff result into "img_full" matrix based on indices of aliased pixels
% % % %                     img_full(ii, aliased_index(index(:,:,:,SlcCount)), SlcCount) = v(StartIndex:EndIndex);
% % % %                     StartIndex = EndIndex+1;
% % % %                     if SlcCount ~= R_slc
% % % %                         EndIndex = EndIndex + sum(index(:,:,:,SlcCount+1));
% % % %                     end
% % % %                 end
% % % %             end
% % % %         end
% % % %     end
 

     
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_1d__general.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
