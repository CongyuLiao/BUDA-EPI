function [img_msl]=recon_SMS_BUDA_CLv2 (kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre_shift, AccY, AccZ,img_applied_topup,img_fieldmap,opt)
PhaseShiftBase= opt.PhaseShiftBase;
phase_diff_flag=opt.phase_diff_flag;  % --0 don't use phase difference --1 use phase differnece for mussels recon
fista =0;
if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
end
num_slc_raw = size(kspace_cor_ap,4);
[N(1), N(2), num_chan, num_slc]= size(sens_gre_shift);
img_msl=zeros(size(img_applied_topup));
    for ii_slc=1:num_slc_raw 
        slc_select = [ii_slc,ii_slc+num_slc/AccZ];
        img_recon_applytopup(:,:,1:2) = cat(3, img_applied_topup(:,:,1,slc_select(1)), img_applied_topup(:,:,2,slc_select(1)));
        img_recon_applytopup(:,:,3:4) = cat(3, img_applied_topup(:,:,1,slc_select(2)), img_applied_topup(:,:,2,slc_select(2)));
        if phase_diff_flag==1
            phase_diff1 = exp(-1i* angle( img_recon_applytopup(:,:,1).*conj(img_recon_applytopup(:,:,2))));
            phase_diff2 = exp(-1i* angle( img_recon_applytopup(:,:,3).*conj(img_recon_applytopup(:,:,2))));
        end
        B0_select1 = 2*pi * img_fieldmap(:,:,slc_select(1));
        B0_select2 = 2*pi * img_fieldmap(:,:,slc_select(2));

         if (PhaseShiftBase ~=0)
            index = 0:AccZ-1 ;
            for jj=1:AccZ
                sensitivity(:,:,:,jj) = circshift(sens_gre_shift(:,:,:,slc_select(jj)),ceil(N(2)/(AccY*peshift))*index(jj),2);
            end
            sens1 = sensitivity(:,:,:,1);
            sens2 = sensitivity(:,:,:,2);
            if phase_diff_flag==1
                phase_diff1= circshift(phase_diff1,ceil(N(2)/(AccY*peshift))*index(1),2);
                phase_diff2= circshift(phase_diff2,ceil(N(2)/(AccY*peshift))*index(2),2);
            end
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
            if phase_diff_flag==1
                EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_pa , cat(2,permute(sens1(:,:,c).* phase_diff1,[3,2,1]), permute(sens2(:,:,c).* phase_diff2,[3 2 1])));    
            else
                EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times,EW_pa , cat(2,permute(sens1(:,:,c),[3,2,1]), permute(sens2(:,:,c),[3 2 1])));
            end
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
 
        clear im_rc
        temp = zeross([N(1), N(2)*2, 2]);
        im_rc = zeross([N(2)*2, 2 , N(1)]);
        
        y_k = im_rc;
        
        x_kn1 = y_k;
        
        t_k = 1;

        for iter = 1:num_iter

%             im_prev = im_rc;

            im_rc = y_k;
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
            mosaic(im_rc_display, 2, 2, 71, num2str(iter), [0,1], 90)

            temp(:,1:end/2,1) = im_rc(:,:,1); % ap slc_select1
            temp(:,1:end/2,2) = im_rc(:,:,2); % pa slc_select1
            temp(:,end/2+1:end,1) = im_rc(:,:,3); % ap slc_select2
            temp(:,end/2+1:end,2) = im_rc(:,:,4); % pa slc_select2
            im_rc = temp;
            im_rc = permute(im_rc, [2,3,1]);

%             update = rmse(im_prev,im_rc)
%             if update < tol
%                 break
%             end
            
            
            
            
            if fista
                t_kp1 = (1 + sqrt(1 + 4 * t_k^2)) / 2;
                coef_kn1 = -(t_k - 1) / t_kp1;
                coef_k = (t_kp1 + t_k - 1) / t_kp1;
            else
                coef_k = 1;
                coef_kn1 = 0;
                t_kp1 = 1;
            end
            x_k= im_rc;
            %update_last = update;
            update = rmse(x_k, x_kn1);
            disp(['iteration: ', num2str(iter), '  update: ', num2str(update), ' %   coef_k: ', num2str(coef_k), '   coef_k-1: ', num2str(coef_kn1)])
            
            y_kp1 = coef_k * x_k + coef_kn1 * x_kn1;
            t_k = t_kp1;
            y_k = y_kp1;
            x_kn1 = x_k;
        
        end
        

        
        
        img_msl(:,:,1,ii_slc) =im_rc_display(:,:,1);
        img_msl(:,:,2,ii_slc) =im_rc_display(:,:,2);
        img_msl(:,:,1,ii_slc+num_slc/AccZ) =im_rc_display(:,:,3);
        img_msl(:,:,2,ii_slc+num_slc/AccZ) =im_rc_display(:,:,4);
    end

end