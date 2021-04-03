function [kspace,ky_idx]=detect_kyLines(K_epi)
 
K_epi= squeeze(K_epi);
% remove readout oversampling
K_epi = ifftc(K_epi,1);
kspace = fftc(K_epi(1+end/4:3*end/4,:,:,:),1);

% detect ky lines
msk = abs(kspace(:,:,1,1)).^(1/3) > 0.001; % Mask to set unacquired points to zero.
ky_idx = find(msk(1,:)~=0);

end