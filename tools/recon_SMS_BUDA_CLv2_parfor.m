function [img_msl]=recon_SMS_BUDA_CLv2_parfor (kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre_shift, AccY, AccZ,img_fieldmap,opt)

%%
PhaseShiftBase= opt.PhaseShiftBase;

if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
else
    peshift = 0;
end
num_slc_raw = size(kspace_cor_ap,4);
[N(1), N(2), num_chan, num_slc]= size(sens_gre_shift);

img_msl=zeros(N(1),N(2),4,num_slc_raw);


%% determine how many  cores u gonna use
 delete(gcp('nocreate'))
 c = parcluster('local'); % build the 'local' cluster object
 total_cores = c.NumWorkers;
 be_mercy=opt.show_mercy;  % how many cores you wanna save for your dear colleagues, show mercy
 
 for ii=1:num_slc_raw
     use_cores=ceil(num_slc_raw/ii);
     if use_cores<=total_cores-be_mercy
         break
     end
 end% parpool(use_cores)%16

 parfor ii_slc=1:num_slc_raw
%for ii_slc = opt.slc
    slc_select = [ii_slc,ii_slc+num_slc/AccZ];
    
    B0_select1 = 2*pi * img_fieldmap(:,:,ii_slc);
    B0_select2 = 2*pi * img_fieldmap(:,:,ii_slc+num_slc/AccZ);
    
    if (PhaseShiftBase ~=0)
        index = 0:AccZ-1 ;
        sens1 = circshift(sens_gre_shift(:,:,:,ii_slc),ceil(N(2)/(AccY*peshift))*index(1),2);
        sens2 = circshift(sens_gre_shift(:,:,:,ii_slc+num_slc/AccZ),ceil(N(2)/(AccY*peshift))*index(2),2);
        
        B0_select1= circshift(B0_select1,ceil(N(2)/(AccY*peshift))*index(1),2);
        B0_select2= circshift(B0_select2,ceil(N(2)/(AccY*peshift))*index(2),2);
    else
        sens1 = sens_gre_shift(:,:,:,slc_select(1));
        sens2 = sens_gre_shift(:,:,:,slc_select(2));
    end
    
    
    
    % AP data
    kspace_ap_collaps = kspace_cor_ap(:,:,:,slc_select(1));
    kspace_ap = kspace_ap_collaps(:,ky_idx_ap,:);
    
    % PA data
    kspace_pa_collaps = kspace_cor_pa(:,:,:,slc_select(1));
    kspace_pa = kspace_pa_collaps(:,ky_idx_pa,:);
    
    signal_ap = ifftc(kspace_ap,1);
    signal_pa = ifftc(kspace_pa,1);
    
    E = fftc(eye(N(2)),1);
    E_ap = E(ky_idx_ap, : );
    E_pa = E(ky_idx_pa, : );
    PE_line = length(ky_idx_ap);
    t_value_ap = [0:AccY:PE_line*AccY-1] * opt.esp;
    t_value_pa = t_value_ap(end:-1:1);
    
    % create and store encoding matrix
    EWC_ap = zeross([num_chan*PE_line, N(2)*2, N(1)]);
    EWC_pa = zeross([num_chan*PE_line, N(2)*2, N(1)]);
    
    
    W_ap1 = exp(1i*mtimesx(repmat(t_value_ap.',[1 1 N(1)]),permute(B0_select1,[3 2 1])));
    W_ap2 = exp(1i*mtimesx(repmat(t_value_ap.',[1 1 N(1)]),permute(B0_select2,[3 2 1])));
    W_pa1 = exp(1i*mtimesx(repmat(t_value_pa.',[1 1 N(1)]),permute(B0_select1,[3 2 1])));
    W_pa2 = exp(1i*mtimesx(repmat(t_value_pa.',[1 1 N(1)]),permute(B0_select2,[3 2 1])));
    EW_ap = bsxfun(@times, repmat(E_ap,1,2),cat(2,W_ap1,W_ap2));
    EW_pa = bsxfun(@times, repmat(E_pa,1,2),cat(2,W_pa1,W_pa2));
    
    for c = 1:num_chan
        EWC_ap(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_ap , cat(2,permute(sens1(:,:,c),[3,2,1]), permute(sens2(:,:,c),[3 2 1])));
        
        EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_pa , cat(2,permute(sens1(:,:,c),[3,2,1]), permute(sens2(:,:,c),[3 2 1])));
        
    end
    
    EWC_apN = mtimesx(EWC_ap,'c',EWC_ap);
    EWC_paN = mtimesx(EWC_pa,'c',EWC_pa);
    rhs_ap = permute(signal_ap,[2 3 1]);
    rhs_ap = reshape(rhs_ap,size(rhs_ap,1)*size(rhs_ap,2),1,size(rhs_ap,3));
    EWC_apHrhs = squeeze(mtimesx(EWC_ap,'c',rhs_ap));
    
    rhs_pa = permute(signal_pa,[2 3 1]);
    rhs_pa = reshape(rhs_pa,size(rhs_pa,1)*size(rhs_pa,2),1,size(rhs_pa,3));
    EWC_paHrhs = squeeze(mtimesx(EWC_pa,'c',rhs_pa));
    
    
    %% mussels recon
    step_size = opt.step_size;
    num_iter = opt.num_iter;
    tol = opt.tol;
    
    winSize = opt.winSize;
    lambda_msl = opt.lambda_msl;
    
    num_sho = 2;
    
    
    temp = zeross([N(1), N(2)*2, 2]);
    im_rc = zeross([N(2)*2, 2 , N(1)]);
    im_rc_sms=zeros(N(1),N(2),4);
    
    for iter = 1:num_iter
        
        im_prev = im_rc;
        
        im_rc(:,1,:) = squeeze(im_rc(:,1,:))- step_size * squeeze(mtimesx(EWC_apN,im_rc(:,1,:)))-EWC_apHrhs;
        im_rc(:,2,:) = squeeze(im_rc(:,2,:))- step_size * squeeze(mtimesx(EWC_paN,im_rc(:,2,:)))-EWC_paHrhs;
        
        
        % mussels constraint
        im_rc = permute(im_rc, [3,1,2]);
        
        im_rc_sms(:,:,1) = im_rc(:,1:end/2,1); % ap slc_select1
        im_rc_sms(:,:,2) = im_rc(:,1:end/2,2); % pa slc_select1
        im_rc_sms(:,:,3) = im_rc(:,end/2+1:end,1); % ap slc_select2
        im_rc_sms(:,:,4) = im_rc(:,end/2+1:end,2); % pa slc_select2
        im_rc = im_rc_sms;
        
        
        % slice 1st
        A = Im2row( fft2call(im_rc(:,:,1:2)), winSize );
        
        [U, S, V] = svdecon(A);
        
        U(isnan(U))=0;S(isnan(S))=0;V(isnan(V))=0;
        U(isinf(U))=0;S(isinf(S))=0;V(isinf(V))=0;
        
        
        keep = 1:floor(lambda_msl*prod(winSize));
        A = U(:,keep) * S(keep,keep) * V(:,keep)';
        
        k_pocs = Row2im(A, [N, num_sho], winSize);
        im_rc(:,:,1:2) = ifft2call(k_pocs);
        
        
        % slice 2nd
        A = Im2row( fft2call(im_rc(:,:,3:4)), winSize );
        
        [U, S, V] = svdecon(A);
        U(isnan(U))=0;S(isnan(S))=0;V(isnan(V))=0;
        U(isinf(U))=0;S(isinf(S))=0;V(isinf(V))=0;
        
        keep = 1:floor(lambda_msl*prod(winSize));
        A = U(:,keep) * S(keep,keep) * V(:,keep)';
        
        k_pocs = Row2im(A, [N, num_sho], winSize);
        im_rc(:,:,3:4) = ifft2call(k_pocs);
        
        %             mosaic(im_rc, 2, 2, 71, num2str(iter), [0,1], 90)
        % shift back
        if (PhaseShiftBase ~=0)
            im_rc_display(:,:,1:2) = circshift(im_rc(:,:,1:2),-ceil(N(2)/(AccY*peshift))*index(1),2);
            im_rc_display(:,:,3:4) = circshift(im_rc(:,:,3:4),-ceil(N(2)/(AccY*peshift))*index(2),2);
        else
            im_rc_display = im_rc;
        end
        %  mosaic(im_rc_display, 2, 2, 71, num2str(iter), [0,1], 90)
        
        temp(:,1:end/2,1) = im_rc(:,:,1); % ap slc_select1
        temp(:,1:end/2,2) = im_rc(:,:,2); % pa slc_select1
        temp(:,end/2+1:end,1) = im_rc(:,:,3); % ap slc_select2
        temp(:,end/2+1:end,2) = im_rc(:,:,4); % pa slc_select2
        im_rc = temp;
        im_rc = permute(im_rc, [2,3,1]);
        
        update = rmse(im_prev,im_rc)
        disp(['iteration: ', num2str(iter), '  update: ', num2str(update), ' % '])
        
        if update < tol
            break
        end
    end
    
    img_msl(:,:,:,ii_slc)=im_rc_display;
end
 delete(gcp('nocreate'))

img_msl=reshape(permute(reshape(img_msl,[N(1),N(2),2,2,num_slc_raw]),[1,2,3,5,4]),[N(1),N(2),2,num_slc]);


end