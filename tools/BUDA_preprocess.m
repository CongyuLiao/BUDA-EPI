function [kspace_cor]= BUDA_preprocess(kspace_cor,prot,ii_rep,kyshift)
pf = prot.ucPhasePartialFourier;
PE_raw= size(kspace_cor,2);

AccY = prot.lAccelFactPE;
PE = ceil(PE_raw/AccY)*AccY;
kspace_cor = mrir_zeropad(kspace_cor,[0 PE-PE_raw 0 0 0 0 0 0 0 0 0],'pre');

% NCO correction for AP/PA shots
 nco_polarity = [-1, +1]; 
 kspace_cor(:,:,:,1,1,1,1,:,1,:) = msEPI_kyshift_correction_v0(kspace_cor(:,:,:,1,1,1,1,:,1,:), prot, kyshift(ii_rep), nco_polarity(ii_rep)) ;


if ii_rep==1 % for ap direction
    if pf < 1
        pad_lines = prot.lPhaseEncodingLines - size(kspace_cor,2);
        kspace_cor = mrir_zeropad(kspace_cor,[0 pad_lines 0 0 0 0 0 0 0 0 0],'pre');
    end
    
    kspace_cor =circshift( kspace_cor ,kyshift(ii_rep),2);
    img = sq(ifft2call(kspace_cor));
    img_crop = img(1+end/4:3*end/4,:,:,:);  % remove readout oversampling
        
else % for pa direction
    kspace_cor = fft2call(flip( ifft2call(kspace_cor), 2));
    if pf < 1
        pad_lines = prot.lPhaseEncodingLines - size(kspace_cor,2);
        kspace_cor = mrir_zeropad(kspace_cor,[0 pad_lines 0 0 0 0 0 0 0 0 0],'post');
    end

    kspace_cor =circshift( kspace_cor ,kyshift(ii_rep),2);
    img = sq(ifft2call(kspace_cor));
    img_crop = img(1+end/4:3*end/4,:,:,:);  % remove readout oversampling
    
end

% mosaic(sq(rsos(kspace_cor(:,:,:,:,:,:,:,:,:,10), 3)), 1, 1, 200+ii_rep, 'k after gc', [0,1]), setGcf(.5)
% mosaic(imrotate(sq(rsos(img_crop, 3)), 90), 3, 5, 210+ii_rep, 'img after gc', [0,1]),



end