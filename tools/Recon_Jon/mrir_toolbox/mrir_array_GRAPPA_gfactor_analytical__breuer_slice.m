function [covmtxFull, g, g_ind, n, n_ind] = mrir_array_GRAPPA_gfactor_analytical__breuer_slice( w, KernalSize_Slice, covmtx, NCol, NLin, varargin)

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
% 3) adopt it for multislice accel case

%%% INPUTS
% w                 grappa kernal information w(:,ch,slc)
% KernalSize_Slice  [kx, ky]
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

% w(:,ch,slc) ---> w_reshaped(kx,ky,ch,ch,slc)
w_reshaped = reshape(w,[KernalSize_Slice(1),KernalSize_Slice(2),nChn,nChn,NslicesEX]);        
w_reshaped = w_reshaped(end:-1:1,end:-1:1,:,:,:);
K = mrir_iDFT(mrir_iDFT(w_reshaped, 1, NCol), 2, NLin);
disp('assume oversample in read by 2')
K = mrir_image_crop(K, 2);

i_noisecov = covmtx * NCol * NLin;
i_noisestd = sqrt(diag(i_noisecov));

covmtxFull = zeros(size(K,1),size(K,2),nChn,nChn,NslicesEX);

for slc = 1:NslicesEX
    slc
    for lin = 1:NLin,
        for col = 1:size(K,1),
            
            % matrix W: [Ctrg x Csrc]
            W = squeeze(K(col, lin, :, :, slc)).';
            
            covmtxFull(col, lin, :, :, slc) = W * i_noisecov * W' ;
            
            if nargout > 1
                n_ind(col, lin, :, slc) = diag(squeeze(sqrt(abs(covmtxFull(col,lin,:,:,slc)))));
                g_ind(col, lin, :, slc) = squeeze(n_ind(col, lin, :, slc)) ./ i_noisestd;
                
                p = squeeze(combine_weights(col,lin,:,1,1,1,1,1,1,slc)).';
                
                n(col, lin, slc) = sqrt(abs( p * W * i_noisecov * W' * p' ));
                g(col, lin, slc) = n(col, lin, slc) ./ sqrt( abs(p * i_noisecov * p') );
            end
            
        end;
    end;
end



return;

