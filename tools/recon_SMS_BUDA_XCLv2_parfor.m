function [img_msl]=recon_SMS_BUDA_XCLv2_parfor (kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre, AccY, AccZ,img_fieldmap,opt)

%%
PhaseShiftBase= opt.PhaseShiftBase;

if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
else
    peshift = 0;
end
num_slc_raw = size(kspace_cor_ap,4);
[N(1), N(2), num_chan, num_slc]= size(sens_gre);
num_k_line=length(ky_idx_ap);

img_msl=zeros(N(1),N(2),2*AccZ,num_slc_raw);


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
end
parpool(use_cores)%16

%parfor
parfor ii_slc=1:num_slc_raw
    slc_select = ii_slc;
    if AccZ>1
        for ii_temp=2:AccZ
            slc_select=[slc_select,ii_slc+num_slc_raw*(ii_temp-1)];
        end
    end
    
    B0_select=zeros(N(1),N(2),AccZ);
    for ii_temp=1:AccZ
        B0_select(:,:,ii_temp)=2*pi * img_fieldmap(:,:,slc_select(ii_temp));
    end
    
    sens=zeros(N(1),N(2),num_chan,AccZ);
    if (PhaseShiftBase ~=0)
        pes_index = mod((1:AccZ)-1,2) ;
        %pes_index=0:AccZ-1;
        % mb2: 0-1
        % mb3: 0-1-0
        % mv4: 0-1-0-1
        for ii_temp=1:AccZ
        sens(:,:,:,ii_temp)=circshift(sens_gre(:,:,:,slc_select(ii_temp)),ceil(N(2)/(AccY*peshift))*pes_index(ii_temp),2);
        B0_select(:,:,ii_temp)= circshift(B0_select(:,:,ii_temp),ceil(N(2)/(AccY*peshift))*pes_index(ii_temp),2);
        end
    else
        for ii_temp=1:AccZ
        sens(:,:,:,ii_temp) = sens_gre(:,:,:,slc_select(ii_temp));
        end
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
    EWC_ap = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
    EWC_pa = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
    
    W_ap=zeross([num_k_line,N(2),N(1),AccZ]);
    W_pa=zeross([num_k_line,N(2),N(1),AccZ]);
    for ii_temp=1:AccZ
    W_ap(:,:,:,ii_temp)=exp(1i*mtimesx(repmat(t_value_ap.',[1 1 N(1)]),permute(B0_select(:,:,ii_temp),[3 2 1])));
    W_pa(:,:,:,ii_temp)=exp(1i*mtimesx(repmat(t_value_pa.',[1 1 N(1)]),permute(B0_select(:,:,ii_temp),[3 2 1])));
    end
    EW_ap = bsxfun(@times, repmat(E_ap,1,AccZ),reshape(permute(W_ap,[1,2,4,3]),[num_k_line,N(2)*AccZ,N(1)]));
    EW_pa = bsxfun(@times, repmat(E_pa,1,AccZ),reshape(permute(W_pa,[1,2,4,3]),[num_k_line,N(2)*AccZ,N(1)]));
     
    for c = 1:num_chan
        EWC_ap(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_ap , reshape(permute(sens(:,:,c,:),[3,2,4,1]),[1,N(2)*AccZ,N(1)]));
        EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_pa , reshape(permute(sens(:,:,c,:),[3,2,4,1]),[1,N(2)*AccZ,N(1)]));         
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
    
    
    im_rc = zeross([N(2)*AccZ, 2 , N(1)]);
  
    for iter = 1:num_iter
        
        im_prev = im_rc;
        
        im_rc(:,1,:) = squeeze(im_rc(:,1,:))- step_size * squeeze(mtimesx(EWC_apN,im_rc(:,1,:)))-EWC_apHrhs;
        im_rc(:,2,:) = squeeze(im_rc(:,2,:))- step_size * squeeze(mtimesx(EWC_paN,im_rc(:,2,:)))-EWC_paHrhs;
        
        
        % mussels constraint
        im_rc = permute(im_rc, [3,1,2]);
             
        im_rc = reshape(permute(reshape(im_rc,[N(1),N(2),AccZ,2]),[1,2,4,3]),[N(1),N(2),AccZ*2]);

        
        for ii_temp=1:AccZ
        % ii slice 
        A = Im2row( fft2call(im_rc(:,:,ii_temp*2-1:ii_temp*2)), winSize );
        
        [U, S, V] = svdecon(A);
        
        U(isnan(U))=0;S(isnan(S))=0;V(isnan(V))=0;
        U(isinf(U))=0;S(isinf(S))=0;V(isinf(V))=0;
        
        
        keep = 1:floor(lambda_msl*prod(winSize));
        A = U(:,keep) * S(keep,keep) * V(:,keep)';
        
        k_pocs = Row2im(A, [N, num_sho], winSize);
        im_rc(:,:,ii_temp*2-1:ii_temp*2) = ifft2call(k_pocs);
        end
        
        
        
        
        %             mosaic(im_rc, 2, 2, 71, num2str(iter), [0,1], 90)
        % shift back
        if (PhaseShiftBase ~=0)
            for ii_temp=1:length(pes_index)
            im_rc_display(:,:,(ii_temp*2-1):ii_temp*2) = circshift(im_rc(:,:,(ii_temp*2-1):ii_temp*2),-ceil(N(2)/(AccY*peshift))*pes_index(ii_temp),2);
            end
         else
            im_rc_display = im_rc;
        end
        %  mosaic(im_rc_display, 2, 2, 71, num2str(iter), [0,1], 90)
        
        im_rc=reshape(permute(reshape(im_rc,[N(1),N(2),2,AccZ]),[1,2,4,3]),[N(1),N(2)*AccZ,2]);
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

img_msl=reshape(permute(reshape(img_msl,[N(1),N(2),2,AccZ,num_slc_raw]),[1,2,3,5,4]),[N(1),N(2),2,num_slc]);


end