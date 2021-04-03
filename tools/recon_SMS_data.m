function [img_recon_pad] = recon_SMS_data (kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre,AccY, AccZ,PhaseShiftBase)
if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
end

[N(1), N(2), num_chan, num_slc]= size(sens_gre);
img_recon_pad = zeross([N,num_slc,2]);
for ii_slc=1:size(kspace_cor_ap,4)
    slc_select = [ii_slc,ii_slc+num_slc/AccZ];
    
    % AP data
    kspace_ap_collaps = kspace_cor_ap(:,:,:,slc_select(1));
    kspace_ap = kspace_ap_collaps(:,ky_idx_ap,:);
    
    % PA data
    kspace_pa_collaps = kspace_cor_pa(:,:,:,slc_select(1));
    kspace_pa = kspace_pa_collaps(:,ky_idx_pa,:);
    
    
    signal_ap = ifftc(kspace_ap,1);
    signal_pa = ifftc(kspace_pa,1);

    
    % sens pre
    sens_gre_shift =sens_gre(:,:,:,slc_select);
    
    if (PhaseShiftBase ~=0)
        sensitivity = zeros(size(sens_gre_shift));
        index = 0:AccZ-1 ;
        for jj=1:AccZ
            sensitivity(:,:,:,jj) = circshift(sens_gre_shift(:,:,:,jj),ceil(N(2)/(AccY*peshift))*index(jj),2);
        end      
        sens1 = sensitivity(:,:,:,1);
        sens2 = sensitivity(:,:,:,2);
    else
        sens1 = sens_gre_shift(:,:,:,1);
        sens2 = sens_gre_shift(:,:,:,2);

    end
    
    
    % sms-sense recon separately w/o B0
    E = fftc(eye(N(2)),1);
    E_ap = E(ky_idx_ap, : );
    E_pa = E(ky_idx_pa, : );
    
    lambda = 1e-3;
    img_sense = zeross([N .* [1,2],2]);
    
    PE_line = length(ky_idx_ap);
    
    for xn = 1:N(1)
        
        svec_ap = signal_ap(xn,:,:);
        svec_ap = svec_ap(:);
        
        svec_pa = signal_pa(xn,:,:);
        svec_pa = svec_pa(:);
        
        EC_ap = zeross([num_chan*PE_line, N(2)*2]);
        EC_pa = zeross([num_chan*PE_line, N(2)*2]);
        
        for c = 1:num_chan
            EC_ap(1 + (c-1)*PE_line : c*PE_line, :) = E_ap * [diag( sens1(xn,:,c)) diag( sens2(xn,:,c))];
            EC_pa(1 + (c-1)*PE_line : c*PE_line, :) = E_pa * [diag( sens1(xn,:,c)) diag( sens2(xn,:,c))];
        end
        
        [U,S,V] = svdecon(EC_ap);
        U(isnan(U))=0;S(isnan(S))=0;V(isnan(V))=0;
        U(isinf(U))=0;S(isinf(S))=0;V(isinf(V))=0;
        img_sense(xn,:,1) = V * diag(diag(S) ./ (diag(S).^2 + lambda)) * U' * svec_ap;
        
        [U,S,V] = svdecon(EC_pa);
        U(isnan(U))=0;S(isnan(S))=0;V(isnan(V))=0;
        U(isinf(U))=0;S(isinf(S))=0;V(isinf(V))=0;
        img_sense(xn,:,2) = V * diag(diag(S) ./ (diag(S).^2 + lambda)) * U' * svec_pa;
    end

        
    img_sense_sms(:,:,1) = img_sense(:,1:end/2,1); % ap select1
    img_sense_sms(:,:,2) = img_sense(:,1:end/2,2); % pa select1
    img_sense_sms(:,:,3) = img_sense(:,end/2+1:end,1); % ap select2
    img_sense_sms(:,:,4) = img_sense(:,end/2+1:end,2); % pa select2

    img_sense = img_sense_sms;

    
    % shift back
    if (PhaseShiftBase ~=0)        
      img_sense(:,:,1:2) = circshift(img_sense(:,:,1:2),-ceil(N(2)/(AccY*peshift))*index(1),2);
      img_sense(:,:,3:4) = circshift(img_sense(:,:,3:4),-ceil(N(2)/(AccY*peshift))*index(2),2);     
    end
    
    %         kspace_sense = fft2call(img_sense);
    mosaic(rot90(img_sense),2,2,12, 'mag ', [0,1])
    mosaic(rot90(angle(img_sense)),2,2,13, 'phase ', [-pi, pi])
    %         mosaic(kspace_sense,2,2,11,'k cplx w/o B0 cor', [0,2]), colormap jet
    
    img_recon_pad(:,:,slc_select(1),:) = img_sense(:,:,1:2); % sense+pf result to estimate phase difference
    img_recon_pad(:,:,slc_select(2),:) = img_sense(:,:,3:4);
    
    clc
    ii_slc
    
end
end