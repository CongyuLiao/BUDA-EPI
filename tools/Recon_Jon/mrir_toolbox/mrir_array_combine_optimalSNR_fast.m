function img = mrir_array_combine_optimalSNR_fast(img_uncombined, sensitivity, covmtx, NormalizationType, covmtxInv)

% This version is very fast because do all pixels in one go instead of
% using 'for' loops

%MRIR_ARRAY_COMBINE_OPTIMALSNR
%
% img_n = mrir_array_combine_optimalSNR(img_uncombined, sensitivity, covmtx, NormalizationType)
%
%
% example:
%
%  covmtx = mrir_array_stats_matrix(double(noise), 'cov');
%  sensitivity = mrir_sensitivity_map(double(raw), 100);
%  img_uncombined = mrir_conventional(double(raw));
%
%  img_n = mrir_array_combine_optimalSNR(img_uncombined, sensitivity, covmtx, 1);

% Kawin Setsompop, August 23 2009
% modified from Jonathan Polimeni code 'mrir_array_combine_opt' 

% Kawin Setsompop, July 9 2010
% modify to make it work for covmtx that is pixel dependent
% covmtx(col, lin, slc, chns, chns)

%**************************************************************************%
Ncol = size(sensitivity, 1);
Nlin = size(sensitivity, 2);
Ncha = size(sensitivity, 3);

% arbitrary constant (normalizes image to approx same level as RSS combo)
%alpha = sqrt(mean(diag(covmtx)));
alpha = 1;

dims = size(img_uncombined);

dims(3) = 1;

% preallocate
img = zeros(dims);

a = mrir_array_combine(sensitivity(:, :, :, 1, 1, 1, 1, 1, 1, :),0); 
mask = double(a>0.5);


if ndims(covmtx) == 2 %i.e. normal covmtx
    
  Psi_inv = inv(covmtx);

  
  for slc = 1:mrir_ice_dimensions(img_uncombined, 'slc'),
      for rep = 1:mrir_ice_dimensions(img_uncombined, 'rep')
          for set = 1:mrir_ice_dimensions(img_uncombined, 'set')
              if rep == 1 && set == 1
                  S = squeeze(sensitivity(:, :, :, 1, 1, 1, 1, 1, 1, slc));
                  S = squeeze(reshape(S,1,[],Ncha)).';
                  Sh = conj(S);
                  W = ((Psi_inv).' * Sh);
                  W_temp = reshape(sum(W.*S,1),Ncol,Nlin);
                  W_temp(W_temp == 0) = 1;
                  if NormalizationType == 1 % not normalize
                      normalization = ones(Ncol,Nlin);
                  elseif NormalizationType == 2 %uniform signal
                      normalization = alpha./ W_temp;
                      %    normalization = alpha./ reshape(sum(W.*S,1),Ncol,Nlin);
                  elseif NormalizationType == 3 %uniform noise
                      normalization = alpha./ sqrt(W_temp);
                      %   normalization = alpha./ sqrt(reshape(sum(W.*S,1),Ncol,Nlin));
                  end
              end
              Ic = squeeze(img_uncombined(:, :, :, set, 1, 1, rep, 1, 1, slc));
              Ic = squeeze(reshape(Ic,1,[],Ncha)).';
              img(:, :, 1, set, 1, 1, rep, 1, 1, slc) = reshape(sum(W.*Ic,1),Ncol,Nlin).*normalization.*mask(:, :, slc);
          end
      end
  end
  
else % covmtx for each voxel location!!
    
    % covmtx(col, lin, chns, chns, slc)
        
    mask = mrir_array_combine(sensitivity,0) > 0;
        
    for slc = 1:mrir_ice_dimensions(img_uncombined, 'slc'),
        for col = 1:size(covmtx,1),
            for lin = 1:size(covmtx,2),
                if mask(col, lin, slc) ~= 0
                    
                    if nargin == 5
                        Psi_inv = squeeze(covmtxInv(col, lin,  :,:, slc));
                    else
                        Psi_inv = inv(squeeze(covmtx(col, lin,  :,:, slc)));
                    end
                    
                    for rep = 1:mrir_ice_dimensions(img_uncombined, 'rep')
                        for set = 1:mrir_ice_dimensions(img_uncombined, 'set')
                            if rep == 1 && set == 1
                                S = squeeze(sensitivity(col, lin, :, 1, 1, 1, 1, 1, 1, slc));
                                %S = squeeze(reshape(S,1,[],Ncha)).';
                                Sh = conj(S);
                                W = ((Psi_inv).' * Sh);
                                W_temp = sum(W.*S,1);
                                %W_temp = reshape(sum(W.*S,1),Ncol,Nlin);
                                %W_temp(W_temp == 0) = 1;
                                if NormalizationType == 1 % not normalize
                                    normalization = 1;
                                elseif NormalizationType == 2 %uniform signal
                                    normalization = alpha/ W_temp;
                                    %    normalization = alpha./ reshape(sum(W.*S,1),Ncol,Nlin);
                                elseif NormalizationType == 3 %uniform noise
                                    normalization = alpha/ sqrt(W_temp);
                                    %   normalization = alpha./ sqrt(reshape(sum(W.*S,1),Ncol,Nlin));
                                end
                            end
                            Ic = squeeze(img_uncombined(col, lin, :, set, 1, 1, rep, 1, 1, slc));
                            %Ic = squeeze(reshape(Ic,1,[],Ncha)).';
                            
                            img(col, lin, 1, set, 1, 1, rep, 1, 1, slc) = sum(W.*Ic,1)*normalization;
                        end
                    end
                    
                else
                    img(col, lin, 1, :, 1, 1, :, 1, 1, slc) = 0;    
                end
                
            end
        end
    end
end










