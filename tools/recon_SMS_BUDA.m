function [img_msl]=recon_SMS_BUDA_CLv1 (kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre_shift, AccY, AccZ,img_applied_topup,img_fieldmap,opt)
PhaseShiftBase= opt.PhaseShiftBase;
phase_diff_flag=0;  % --0 don't use phase difference --1 use phase differnece for mussels recon

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
         %   sensitivity = zeros(size(sens_gre_shift));
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

        EWC_apH = zeross([N(2)*2, num_chan*PE_line, N(1)]);       
        EWC_paH = zeross([N(2)*2, num_chan*PE_line, N(1)]);
        
        
        EWC_apN = zeros(size(EWC_ap, 2), size(EWC_ap, 2), size(EWC_ap, 3));
        EWC_paN = zeros(size(EWC_pa, 2), size(EWC_pa, 2), size(EWC_pa, 3));
        
        EWC_apHrhs =  zeros(size(EWC_apH, 1), size(EWC_apH, 3));
        EWC_paHrhs =  zeros(size(EWC_paH, 1), size(EWC_paH, 3));


        for xn = 1:N(1)

            W_ap1 = exp(1i * t_value_ap.' * B0_select1(xn,:));
            W_ap2 = exp(1i * t_value_ap.' * B0_select2(xn,:));
            W_pa1 = exp(1i * t_value_pa.' * B0_select1(xn,:));
            W_pa2 = exp(1i * t_value_pa.' * B0_select2(xn,:));

            EW_ap = repmat(E_ap,1,2) .* [W_ap1 W_ap2];
            EW_pa = repmat(E_pa,1,2) .* [W_pa1 W_pa2];

            for c = 1:num_chan
                EWC_ap(1 + (c-1)*PE_line : c*PE_line, :, xn) = EW_ap * [diag( sens1(xn,:,c)) zeros(N); zeros(N) diag( sens2(xn,:,c))];
                if phase_diff_flag==1
                    EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, xn) = EW_pa * [diag( sens1(xn,:,c).* phase_diff1(xn,:)) zeros(N); zeros(N) diag( sens2(xn,:,c).* phase_diff2(xn,:))];
                else
                    EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, xn) = EW_pa * [diag( sens1(xn,:,c)) zeros(N); zeros(N) diag( sens2(xn,:,c))];
                end
            end

            EWC_apH(:,:,xn) = EWC_ap(:,:,xn)';
            EWC_paH(:,:,xn) = EWC_pa(:,:,xn)';
            EWC_apN(:, :, xn) = EWC_apH(:,:,xn) * EWC_ap(:,:,xn);
            EWC_paN(:, :, xn) = EWC_paH(:,:,xn) * EWC_pa(:,:,xn);
            

            rhs_ap = signal_ap(xn,:,:);
            EWC_apHrhs(:,xn) = EWC_apH(:,:,xn) * rhs_ap(:);
            
            rhs_pa = signal_pa(xn,:,:);
            EWC_paHrhs(:,xn) = EWC_paH(:,:,xn) * rhs_pa(:);
        end

%% mussels recon
        step_size = opt.step_size;
        num_iter = opt.num_iter;
        tol = opt.tol;

        winSize = opt.winSize;
        lambda_msl = opt.lambda_msl;

        num_sho = 2;


        clear im_rc
        im_rc(:,:,1)=cat(2,img_recon_applytopup(:,:,1),img_recon_applytopup(:,:,3));
        im_rc(:,:,2)=cat(2,img_recon_applytopup(:,:,2),img_recon_applytopup(:,:,4));
        im_rc = permute(im_rc, [2,3,1]);
        temp = zeross([N(1), N(2)*2, 2]);
        im_rc = zeross(size(im_rc));

        for iter = 1:num_iter
            im_prev = im_rc;

            for xn = 1:N(1)
                im_rc(:,1,xn) = im_rc(:,1,xn) - step_size * ( EWC_apN(:,:,xn) * im_rc(:,1,xn) - EWC_apHrhs(:, xn) );
                im_rc(:,2,xn) = im_rc(:,2,xn) - step_size * ( EWC_paN(:,:,xn) * im_rc(:,2,xn) - EWC_paHrhs(:, xn) );
            end
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



            update = rmse(im_prev,im_rc)
            if update < tol
                break
            end
        end


        img_msl(:,:,1,ii_slc) =im_rc_display(:,:,1);
        img_msl(:,:,2,ii_slc) =im_rc_display(:,:,2);
        img_msl(:,:,1,ii_slc+num_slc/AccZ) =im_rc_display(:,:,3);
        img_msl(:,:,2,ii_slc+num_slc/AccZ) =im_rc_display(:,:,4);
    end

end