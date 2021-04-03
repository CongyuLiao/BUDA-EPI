function [covmtxFull, g, g_ind, n, n_ind, InvCovmtxFull] = mrir_array_GRAPPA_gfactor_analytical__breuer_Full_v2( w, KernalSize_Slice, Phase, G, R_inplane, covmtx, NCol, NLin, varargin)

% Based on MRIR_ARRAY_GRAPPA_GFACTOR_ANALYTICAL__BREUER written by jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/13
% using algorithm from Breuer et al. (2008), "A general formulation for quantitative g-factor
% calculation in GRAPPA reconstructions [Abstract]"; Proc. ISMRM 16:10.

% kawin setsompop
% July 14 2010
% 1) fix bug in noise calculation
% 2) sort out scaling issue so work with P.F. and account for reduced noise due to inplane acceleration
%  -> explaination:  
%   2.1) noise is scaled when convert to image domain as jon's code uses a non-unity Fourier transform
%   2.2) P.F. doesnt change the scaling of the noise when convert to image domain as in jon's code for P.F. recon, he uses a
%   intensity correction factor to scale up the signal (and hence the noise) of P.F. data to equivilent full data case
% 3) do for both multislice and inplane accel here 

%%% INPUTS
% w                 slice-grappa kernal information w(:,ch,slc)
% KernalSize_Slice  [kx, ky]
% Phase             blip phase shift vector
% G                 'array' of inplane-grappa kernal information (struct)
% R_inplane         inplane accel factor
% covmtx            original coil noise covariance
% NCol              number of col points of 'final' k-space after P.F. and IPAT recon 
% NLin              number of lines points of 'final' k-space after P.F. and IPAT recon 
% combine_weights   this is normally conj(sensitivity profile)
% |->(varargin) 
%
%%% OUTPUTS
% covmtxFull        covmtx for every voxel in the output image (Col,Lin,Ch,Ch,Slc)
% g                 g-factor final recon image
% g_ind             g-factor each coil
% n                 noise 
% n_ind             noise for each coil image

%**************************************************************************%
if length(varargin)   == 1
    combine_weights = varargin{1};
end

if nargout > 1
    if isempty(varargin)
        disp('need to supply combine weights to be able to calc g factors')
        keyboard
    end
end

%**************************************************************************%
NslicesEX = size(w,3);
nChn = length(covmtx);

i_noisecov = covmtx * NCol * NLin;
i_noisecov_accel = i_noisecov/R_inplane;
i_noisestd = sqrt(diag(i_noisecov));

%% slice-Grappa

% w(:,ch,slc) ---> w_reshaped(kx,ky,ch,ch,slc)
w_reshaped = reshape(w,KernalSize_Slice(1),KernalSize_Slice(2),nChn,nChn,NslicesEX);        
w_reshaped = w_reshaped(end:-1:1,end:-1:1,:,:,:);
s = size(w_reshaped);
w_reshaped_full = zeros(s(1),s(2)+ (s(2)-1)*(R_inplane-1),s(3),s(4),s(5));
w_reshaped_full(:,1:R_inplane:end,:,:,:) = w_reshaped;
K = mrir_iDFT(mrir_iDFT(w_reshaped_full, 1, NCol), 2, NLin);
disp('assume oversample in read by 2')
K = mrir_image_crop(K, 2);

covmtxFull_step1 = zeros(size(K,1),size(K,2),nChn,nChn,NslicesEX);

for slc = 1:NslicesEX
    for lin = 1:NLin,
        for col = 1:size(K,1),            
            W = squeeze(K(col, lin, :, :, slc)).';            
            covmtxFull_step1(col, lin, :, :, slc) = W * i_noisecov_accel * W' ;            
        end
    end
end

for SlcCount = 1:NslicesEX % shift back the covmtx so final image is unshifted. 
    covmtxFull_step1(:,:,:,:,SlcCount) = mrir_iDFT_freqencode( mrir_iDFT_phasencode (CaipirinhaShift_K(mrir_fDFT_freqencode( mrir_fDFT_phasencode (covmtxFull_step1(:,:,:,:,SlcCount) )), -Phase(SlcCount)) ));
end


%% inplane accel
covmtxFull = zeros(size(covmtxFull_step1));

combineMethod = 1;
%Ncut = 1;
CondFac = 2 % the larger this is the less we are going to truncate (allow higher cond number)
i_noisecov_INV = inv(i_noisecov);

if combineMethod == 3
    InvCovmtxFull = zeros(size(covmtxFull));
    cond_covmtx = cond(i_noisecov);
else 
    InvCovmtxFull = [];
end

for slc = 1:NslicesEX
    
    k = mrir_array_GRAPPA_conv_kernel(G{slc}); 
    K = mrir_iDFT(mrir_iDFT(k, 1, NCol), 2, NLin);
    disp('assume oversample in read by 2')
    K = mrir_image_crop(K, 2);
    
    for lin = 1:NLin,
        for col = 1:size(K,1),
            W = squeeze(K(col, lin, 1, :, :)).';
            covmtxFull(col, lin, :,:,slc) = W * squeeze(covmtxFull_step1(col, lin, :, :, slc)) * W' ;
            if nargout > 1
                n_ind(col, lin, :, slc) = diag(squeeze(sqrt(abs(covmtxFull(col,lin,:,:,slc)))));
                g_ind(col, lin, :, slc) = squeeze(n_ind(col, lin, :, slc)) ./ i_noisestd /sqrt(R_inplane);
                
                if combineMethod == 1
                    p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1,slc)).';
                    n(col, lin, slc) = sqrt(abs( p * squeeze(covmtxFull(col, lin, :,:,slc)) * p' ));
                    g(col, lin, slc) = n(col, lin, slc) ./ sqrt( abs(p * i_noisecov * p') ) /sqrt(R_inplane);
                elseif combineMethod == 2 % old covariance
                    p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1,slc)).' * i_noisecov_INV;
                    p = p/(p*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1,slc)))); % normalization
                    n(col, lin, slc) = sqrt(abs( p * squeeze(covmtxFull(col, lin, :,:,slc)) * p' ));
                    g(col, lin, slc) = n(col, lin, slc) ./ sqrt( abs(p * i_noisecov * p') ) /sqrt(R_inplane);
                elseif combineMethod == 3 % using covmtx_FUll
                    [U,S,V] = svd(squeeze(covmtxFull(col,lin,:,:,slc)));
                    t = diag(1./diag(S));
                    %t(end-(Ncut-1):end,:) = 0;
                    t(t> [(cond_covmtx*CondFac)*t(1,1)]) = 0;
                    InvCovmtxFull(col,lin,:,:,slc) = V*t*U';
                    p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1,slc)).' * squeeze(InvCovmtxFull(col, lin, :,:,slc));
                    p2 = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1,slc)).' * i_noisecov_INV;
                    p = p/(p*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1,slc)))); % normalization so uniform and hence same signal in both cases
                    p2 = p2/(p2*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1,slc)))); % normalization
                    n(col, lin,slc) = sqrt(abs( p * squeeze(covmtxFull(col, lin, :,:,slc)) * p' ));
                    g(col, lin,slc) = n(col, lin) ./ sqrt( abs(p2 * i_noisecov * p2') ) /sqrt(R_inplane);
                end
                
            end
        end
    end

end


  

