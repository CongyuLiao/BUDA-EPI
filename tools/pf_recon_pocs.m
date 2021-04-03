function [img_recon,img_recon_complex] = pf_recon_pocs(kspace_pk,Nsize,flag)

PE_pk = Nsize(1);       % partial kspace size along PE
PE_full = Nsize(2);     % full kspace size along PE
FE = Nsize(3);          % size along FE (readout)

sym_acq_region  = (PE_full - PE_pk + 1):PE_pk;

if flag == 1
    acq_region      = (PE_full - PE_pk + 1):PE_full;  % for ap pf data
else
    acq_region      = 1:PE_pk;                        % for pa pf data
end

% zero padding
kspace_padding = zeros(FE, PE_full);
kspace_padding(:,acq_region) = kspace_pk;


% phase estimation from symmetric part of data
hamming_filterPE  = zeros(FE, PE_full);
PE_filterLength = length(sym_acq_region);
PE_filterIndex = [-floor(PE_filterLength/2)+1:floor(PE_filterLength/2)] + round(PE_full/2);

hamming_filterPE(:, PE_filterIndex) = repmat(hamming(numel(PE_filterIndex)).', [FE, 1]);
 
 
hamming_filterFE = zeros(FE, PE_full);
RO_filterLength = length(sym_acq_region)*(FE/PE_full);
RO_filterIndex = [-floor(RO_filterLength/2)+1:floor(RO_filterLength/2)] + round(FE/2);
 
hamming_filterFE(RO_filterIndex,:) = repmat(hamming(numel(RO_filterIndex)), [1, PE_full]);
 

hamming_filter =  hamming_filterPE.* hamming_filterFE;


image_filter = ifft2c(hamming_filter.*kspace_padding); 
phase_filter = angle(image_filter);  % estimating image phase 

% phase corrected partial image data
img_padding = ifft2c(kspace_padding); 
img_recon = img_padding .* exp(-1i*phase_filter);



% pocs recon
pocs_tol = 1e-3;               % tolerance

for k=1:100
    img_recon_last = img_recon;
    
    % project onto the data constraint set
    kspace_recon = fft2c(img_recon);

    kspace_recon(:,acq_region) = kspace_pk;

    img_recon = ifft2c(kspace_recon);

    % project onto the phase constraint set
    img_recon = abs(img_recon).*exp(1i*phase_filter);
    
    rela_change = norm2(img_recon_last - img_recon) / norm2(img_recon);
    if rela_change < pocs_tol
        break
    end
    
end

img_recon_complex = img_recon;

% kspace_recon_complex = fft2c(img_recon_complex);
% mosaic(kspace_recon_complex,1,1,1,['pocs iter: ', num2str(k) '  relative change: ', num2str(rela_change)], [0,1e-3]), set(gcf,'color','w')
% mosaic(rot90(img_recon_complex_crop),1,1,2, 'pocs', [0,1.5e-3])


if (1) % do real diffussion
    kspace_recon = fft2c(img_recon);
    kspace_recon(:, acq_region) = kspace_pk;

    img_recon = ifft2c(kspace_recon);
    img_recon = img_recon.*exp(-1i*phase_filter);
end

% kspace_recon = fft2c(img_recon);
% mosaic(kspace_recon,1,1,3,'real k', [0,1e-3]), set(gcf,'color','w')
% mosaic(rot90(img_recon_crop),1,1,4, 'real', [0,1.5e-3])

end






