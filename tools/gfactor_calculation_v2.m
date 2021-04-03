function [gfactor]=gfactor_calculation_v2 (ky_idx_ap,ky_idx_pa,sens_gre_shift, AccY, AccZ,img_fieldmap,opt)

%%
PhaseShiftBase= opt.PhaseShiftBase;

if (PhaseShiftBase(1) ~=0  || PhaseShiftBase(2) ~=0)
    peshift = ceil(2*pi./PhaseShiftBase);
    peshift(isinf(peshift))=0;
else
    peshift = [0,0];
end

[N(1), N(2), num_chan, num_slc]= size(sens_gre_shift);
num_slc_raw = num_slc/AccZ;
% img_msl=zeros(N(1),N(2),4,num_slc_raw);
gfactor=zeros(N(1),N(2)*AccZ,num_slc_raw);

% gfactor=zeros(360,720,21);
%% determine how many  cores u gonna use
% delete(gcp('nocreate'))
% c = parcluster('local'); % build the 'local' cluster object
% total_cores = c.NumWorkers;
% be_mercy=opt.show_mercy;  % how many cores you wanna save for your dear colleagues, show mercy
% 
% for ii=1:num_slc_raw
%     use_cores=ceil(num_slc_raw/ii);
%     if use_cores<=total_cores-be_mercy
%         break
%     end
% end
% parpool(use_cores)%16

for ii_slc=opt.nslc %1:num_slc_raw
    ii_slc
    if (AccZ>1)
        slc_select = [ii_slc,ii_slc+num_slc/AccZ];
        
        B0_select1 = 2*pi * img_fieldmap(:,:,ii_slc);
        B0_select2 = 2*pi * img_fieldmap(:,:,ii_slc+num_slc/AccZ);
        
        if (peshift(1) ~=0  || peshift(2) ~=0)
            index = 0:AccZ-1 ;
            if (peshift(1)~=0)
                sens1 = circshift(sens_gre_shift(:,:,:,ii_slc),ceil(N(2)/(AccY*peshift(1)))*index(1),2);
                B0_select1= circshift(B0_select1,ceil(N(2)/(AccY*peshift(1)))*index(1),2);
            else
                sens1 = sens_gre_shift(:,:,:,slc_select(1));
            end
            if (peshift(2)~=0)
                sens2 = circshift(sens_gre_shift(:,:,:,ii_slc+num_slc/AccZ),ceil(N(2)/(AccY*peshift(2)))*index(2),2);
                B0_select2= circshift(B0_select2,ceil(N(2)/(AccY*peshift(2)))*index(2),2);
            else
                sens2 = sens_gre_shift(:,:,:,slc_select(2));
            end

        else
            sens1 = sens_gre_shift(:,:,:,slc_select(1));
            sens2 = sens_gre_shift(:,:,:,slc_select(2));
        end
        
        E = fftc(eye(N(2)),1);
        E_ap = E(ky_idx_ap, : );
        E_pa = E(ky_idx_pa, : );
        PE_line = length(ky_idx_ap);
        t_value_ap = [0:AccY:PE_line*AccY-1] * opt.esp;
        t_value_pa = t_value_ap(end:-1:1);
        
        % create and store encoding matrix
        EWC_ap = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
        EWC_pa = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
        
        
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
        
        Enc = cat(1, EWC_ap, EWC_pa);
        
        % gfactor estimation by SID
        scs = mtimesx(Enc,'c',Enc);
        tmp = zeros (size(scs,3),size(scs,1));
        for xn = 1: size(scs,3)
            scsx = scs(:,:,xn);
            scsxi =mypinv(scs(:,:,xn));
            for yy = 1:size(scs,1)
                tmp (xn,yy) = abs(sqrt(scsxi(yy, yy) * scsx(yy, yy)));
            end
        end
        gfactor(:,:,ii_slc) = tmp;
%         mosaic(1./tmp, 1, 1, 514, ['1/gfactor map: slice group', num2str(ii_slc)], [0,1], 90),  colormap jet
        
    else
        
        
        B0_select1 = 2*pi * img_fieldmap(:,:,ii_slc);
        sens1 = sens_gre_shift(:,:,:,ii_slc);
   
        E = fftc(eye(N(2)),1);
        E_ap = E(ky_idx_ap, : );
        E_pa = E(ky_idx_pa, : );
        PE_line = length(ky_idx_ap);
        t_value_ap = [0:AccY:PE_line*AccY-1] * opt.esp;
        t_value_pa = t_value_ap(end:-1:1);
        
        % create and store encoding matrix
        EWC_ap = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
        EWC_pa = zeross([num_chan*PE_line, N(2)*AccZ, N(1)]);
        
        
        W_ap1 = exp(1i*mtimesx(repmat(t_value_ap.',[1 1 N(1)]),permute(B0_select1,[3 2 1])));
        W_pa1 = exp(1i*mtimesx(repmat(t_value_pa.',[1 1 N(1)]),permute(B0_select1,[3 2 1])));
        EW_ap = bsxfun(@times, E_ap,W_ap1);
        EW_pa = bsxfun(@times, E_pa,W_pa1);
        
        for c = 1:num_chan
            EWC_ap(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times, EW_ap , permute(sens1(:,:,c),[3,2,1]));          
            EWC_pa(1 + (c-1)*PE_line : c*PE_line, :, :) = bsxfun(@times, EW_pa , permute(sens1(:,:,c),[3,2,1]));
            
        end
     
        Enc = cat(1, EWC_ap, EWC_pa);
        
        % gfactor estimation by SID
        scs = mtimesx(Enc,'c',Enc);
        tmp = zeros (size(scs,3),size(scs,1));
        for xn = 1: size(scs,3)
            scsx = scs(:,:,xn);
            scsxi =mypinv(scs(:,:,xn));
            for yy = 1:size(scs,1)
                tmp (xn,yy) = abs(sqrt(scsxi(yy, yy) * scsx(yy, yy)));
            end
        end
        gfactor(:,:,ii_slc) = tmp;
%         mosaic(1./tmp, 1, 1, 514, ['1/gfactor map: slice group', num2str(ii_slc)], [0,1], 90),  colormap jet
        
        toc
    
    end

end
    % delete(gcp('nocreate'))
end