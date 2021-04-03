function [img_recon_pad] = recon_SMS_data_CLv1(kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre,AccY, AccZ,PhaseShiftBase)
if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
end

[N(1), N(2), num_chan, num_slc]= size(sens_gre);
img_recon_pad = zeross([N,num_slc,2]);

lsqr_iter = 300;
lsqr_tol = 1e-3;
mask_ap = zeros(N(1),N(2),num_chan);
mask_pa = zeros(N(1),N(2),num_chan);
mask_ap(:,ky_idx_ap,:) =1;
mask_pa(:,ky_idx_pa,:) =1;
for ii_slc=1:size(kspace_cor_ap,4)
    slc_select = [ii_slc,ii_slc+num_slc/AccZ];
    
    img_sense = zeros([N.*[AccZ,1],2]);
    
    % sens shift
    sens_gre_shift =sens_gre(:,:,:,slc_select);
    sensitivity = zeros(size(sens_gre_shift));
    if (PhaseShiftBase ~=0)       
        index = 0:AccZ-1 ;
        for jj=1:AccZ
            sensitivity(:,:,:,jj) = circshift(sens_gre_shift(:,:,:,jj),ceil(N(2)/(AccY*peshift))*index(jj),2);
        end      
    else
        sensitivity  = sens_gre_shift;
    end
    
    % AP data
    kspace_ap_collaps = kspace_cor_ap(:,:,:,slc_select(1));
    img_ap_collaps = ifft2call (kspace_ap_collaps );
    img_ap_collaps = repmat(img_ap_collaps, [AccZ,1,1,1,1])/AccZ;
    kspace_ap = fft2call (img_ap_collaps);
    kspace_ap(2:AccZ:end,:,:,:,:) = 0; 
    kspace_ap = kspace_ap.* repmat(mask_ap,[AccZ,1,1,1,1]);
    
    % PA data
    kspace_pa_collaps = kspace_cor_pa(:,:,:,slc_select(1));
    img_pa_collaps = ifft2call (kspace_pa_collaps );
    img_pa_collaps = repmat(img_pa_collaps, [AccZ,1,1,1,1]) /AccZ;
    kspace_pa = fft2call (img_pa_collaps);
    kspace_pa(2:AccZ:end,:,:,:,:) = 0; 
    kspace_pa = kspace_pa.* repmat(mask_pa,[AccZ,1,1,1,1]);
    
    % SMS-sense using LSQR
    param = [];
    param.N = N .* [AccZ,1];
    param.num_chan = num_chan;
    param.lambda = 1e-3;
    param.sens = permute(reshape(permute(sensitivity,[2 3 1 4]),[N(2), num_chan, N(1)*AccZ]),[3 1 2]);
    param.m2d = kspace_ap~=0;  
    % sms-sense for AP shot
    res_ap = lsqr(@apply_sense_tikc, cat(1, kspace_ap(:), zeross([prod(N)*AccZ,1])), lsqr_tol, lsqr_iter, [], [], [], param);
    img_sense(:,:,1) = reshape(res_ap, N.*[AccZ,1]);
    
    param.m2d = kspace_pa~=0;
    % sms-sense for PA shot
    res_pa = lsqr(@apply_sense_tikc, cat(1, kspace_pa(:), zeross([prod(N)*AccZ,1])), lsqr_tol, lsqr_iter, [], [], [], param);
    img_sense(:,:,2) = reshape(res_pa, N.*[AccZ,1]);
    
      
    img_sense_sms(:,:,1) = img_sense(1:end/2,:,1); % ap select1
    img_sense_sms(:,:,2) = img_sense(1:end/2,:,2); % pa select1
    img_sense_sms(:,:,3) = img_sense(end/2+1:end,:,1); % ap select2
    img_sense_sms(:,:,4) = img_sense(end/2+1:end,:,2); % pa select2
    img_sense = img_sense_sms;

    
    % shift back
    if (PhaseShiftBase ~=0)        
      img_sense(:,:,1:2) = circshift(img_sense(:,:,1:2),-ceil(N(2)/(AccY*peshift))*index(1),2);
      img_sense(:,:,3:4) = circshift(img_sense(:,:,3:4),-ceil(N(2)/(AccY*peshift))*index(2),2);     
    end
    
  
%     mosaic(rot90(img_sense),2,2,12, 'mag ', [0,1])
%     mosaic(rot90(angle(img_sense)),2,2,13, 'phase ', [-pi, pi])
    
    img_recon_pad(:,:,slc_select(1),:) = img_sense(:,:,1:2); 
    img_recon_pad(:,:,slc_select(2),:) = img_sense(:,:,3:4);
    
   
    
end
end