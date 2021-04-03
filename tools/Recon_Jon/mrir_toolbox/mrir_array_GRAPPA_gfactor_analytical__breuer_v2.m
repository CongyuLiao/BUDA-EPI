function [covmtxFull, g, g_ind, n, n_ind, InvCovmtxFull] = mrir_array_GRAPPA_gfactor_analytical__breuer_v2(G, covmtx, R_inplane, NCol, NLin, varargin)

% Based on MRIR_ARRAY_GRAPPA_GFACTOR_ANALYTICAL__BREUER written by jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/13
% using algorithm from Breuer et al. (2008), "A general formulation for quantitative g-factor
% calculation in GRAPPA reconstructions [Abstract]"; Proc. ISMRM 16:10.

% kawin setsompop
% July 13 2010
% 1) fix bug in noise calculation
% 2) sort out scaling issue so work with P.F. and account for reduced noise due to inplane acceleration
%  -> explaination:  
%   2.1) noise is scaled when convert to image domain as jon's code uses a non-unity Fourier transform
%   2.2) P.F. doesnt change the scaling of the noise when convert to image domain as in jon's code for P.F. recon, he uses a
%   intensity correction factor to scale up the signal (and hence the noise) of P.F. data to equivilent full data case
%

%%% INPUTS
% G                 grappa kernal information
% covmtx            original coil noise covariance
% NCol              number of col points of 'final' k-space after P.F. and IPAT recon 
% NLin              number of lines points of 'final' k-space after P.F. and IPAT recon 
% combine_weights   this is normally conj(sensitivity profile)
% |->(varargin) 
%
%%% OUTPUTS
% covmtxFull        covmtx for every voxel in the output image (Col,Lin,Ch,Ch)
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

nCha = unique(size(covmtx));

k = mrir_array_GRAPPA_conv_kernel(G);
K = mrir_iDFT(mrir_iDFT(k, 1, NCol), 2, NLin);
disp('assume oversample in read by 2')
K = mrir_image_crop(K, 2);

i_noisecov = covmtx * NCol * NLin;
i_noisecov_IPAT = i_noisecov/R_inplane;
i_noisestd = sqrt(diag(i_noisecov));

covmtxFull = zeros(size(K,1),size(K,2),nCha,nCha);

combineMethod = 3; 
%Ncut = 1;
CondFac = 2 % the larger this is the less we are going to truncate (allow higher cond number)
i_noisecov_INV = inv(i_noisecov);
if combineMethod == 3
    InvCovmtxFull = zeros(size(covmtxFull));
    cond_covmtx = cond(i_noisecov);
else
    InvCovmtxFull = [];
end

for lin = 1:NLin,
    for col = 1:size(K,1),
        
        % matrix W: [Ctrg x Csrc]
        W = squeeze(K(col, lin, 1, :, :)).';
                
        covmtxFull(col, lin, :,:) = W * i_noisecov_IPAT * W' ;
        
        if nargout > 1
            
            n_ind(col, lin, :) = diag(squeeze(sqrt(abs(covmtxFull(col,lin,:,:)))));
            g_ind(col, lin, :) = squeeze(n_ind(col, lin, :)) ./ i_noisestd /sqrt(R_inplane);
            
            if combineMethod == 1 % no covariance
                % vector p: [1 x Ctrg]        1    2  3  4 5 6 7 8    9
                p = squeeze(combine_weights(col, lin, :, 1,1,1,1,1,1)).';
                %n(col, lin) = sqrt(abs( p * W * i_noisecov_IPAT * W' * p' ));
                n(col, lin) = sqrt(abs( p * squeeze(covmtxFull(col,lin,:,:)) * p' ));
                g(col, lin) = n(col, lin) ./ sqrt( abs(p * i_noisecov * p') ) /sqrt(R_inplane);
            elseif combineMethod == 2 % old covariance
                p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1)).' * i_noisecov_INV;
                p = p/(p*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1)))); % normalization
                n(col, lin) = sqrt(abs( p * squeeze(covmtxFull(col, lin, :,:)) * p' ));
                g(col, lin) = n(col, lin) ./ sqrt( abs(p * i_noisecov * p') ) /sqrt(R_inplane);
            elseif combineMethod == 3 % using covmtx_FUll
                [U,S,V] = svd(squeeze(covmtxFull(col,lin,:,:)));
                t = diag(1./diag(S));
                %t(end-(Ncut-1):end,:) = 0;
                t(t> [(cond_covmtx*CondFac)*t(1,1)]) = 0;
                InvCovmtxFull(col,lin,:,:) = V*t*U';               
                p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1)).' * squeeze(InvCovmtxFull(col, lin, :,:));
                p2 = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1)).' * i_noisecov_INV;
                p = p/(p*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1)))); % normalization so uniform and hence same signal in both cases
                p2 = p2/(p2*squeeze(conj(combine_weights(col,lin,:,1,1,1,1,1,1)))); % normalization
                n(col, lin) = sqrt(abs( p * squeeze(covmtxFull(col, lin, :,:)) * p' ));
                g(col, lin) = n(col, lin) ./ sqrt( abs(p2 * i_noisecov * p2') ) /sqrt(R_inplane);                
            end
        end
        
    end;
end;



return;

