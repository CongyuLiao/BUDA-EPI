function [I_training_slice_hybridPh, K_training_slice_hybridPh]=estimate_background_phase_fleet (Img_EPI_comb,Img_ref_comb,Img_ref,SensMap,PhaseDiff_LPF,Irecon_SMS_bkgEst)
%% 2017/09/14 by Congyu Liao
% use complex division to get phase differences 
% used for FLEET ref
%%

% a_setpath_v2;
nCoil=size(Img_ref,3);
AccZ=size(Img_ref,4);


I_PhaseEst_LPF = permute(abs(Img_EPI_comb),[1 2 4 3]).*exp(1i*PhaseDiff_LPF); % use abs to help weight LPF

I_training_slice_hybridPh = abs(Img_ref).*exp(1i*angle(Irecon_SMS_bkgEst));
K_training_slice_hybridPh = mrir_fDFT_freqencode(mrir_fDFT_phasencode(I_training_slice_hybridPh));

for ss = 1:AccZ
    sens_slice = SensMap(:,:,:,ss);
    img_slice = I_training_slice_hybridPh(:,:,:,ss);
    img_combo = sum(conj(sens_slice) .* img_slice,3) ./ (eps + sum(abs(sens_slice).^2,3));
    I_training_slice_hybridPh_comb(:,:,ss) = img_combo;
    bkgPhase(:,:,ss)= angle(I_training_slice_hybridPh_comb(:,:,ss)); %% bkgphase is the difference between SE-dw EPI and SE acs data
end

%% hard coded it
if(0)
figure;axis equal;
if AccZ==1
    % only one slice
%     subplot(1,5,1);imagesc(angle(I_PhaseEst(:,:,1,1))); title('EPI-Ref'); %axis equal;square; %colormap gray;
    subplot(1,5,2);imagesc(angle(Img_ref_comb(:,:,1))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(1,5,3);imagesc(angle(Img_EPI_comb(:,:,1))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(1,5,4);imagesc(angle(I_PhaseEst_LPF(:,:,1,1))); title('EPI-Ref after LPF'); %axis equal;square; %colormap gray;
    subplot(1,5,5);imagesc(angle(I_training_slice_hybridPh_comb(:,:,1))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
elseif AccZ==2
    % Slice 1
%     subplot(2,5,1);imagesc(angle(I_PhaseEst(:,:,1,1))); title('Phase Est'); %axis equal;square; %colormap gray;
    subplot(2,5,2);imagesc(angle(Img_ref_comb(:,:,1))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(2,5,3);imagesc(angle(Img_EPI_comb(:,:,1))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(2,5,4);imagesc(angle(I_PhaseEst_LPF(:,:,1,1))); title('Phase Est after LPF'); %axis equal;square; %colormap gray;
    subplot(2,5,5);imagesc(angle(I_training_slice_hybridPh_comb(:,:,1))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
    % Slice 2 
%     subplot(2,5,6);imagesc(angle(I_PhaseEst(:,:,1,2))); title('Phase Est'); %axis equal;square; %colormap gray;
    subplot(2,5,7);imagesc(angle(Img_ref_comb(:,:,2))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(2,5,8);imagesc(angle(Img_EPI_comb(:,:,2))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(2,5,9);imagesc(angle(I_PhaseEst_LPF(:,:,1,2))); title('Phase Est after LPF'); %axis equal;square; %colormap gray;
    subplot(2,5,10);imagesc(angle(I_training_slice_hybridPh_comb(:,:,2))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
   
elseif AccZ==3
    % Slice 1
%     subplot(3,5,1);imagesc(angle(I_PhaseEst(:,:,1,1))); title('Phase Est'); %axis equal;square; %colormap gray;
    subplot(3,5,2);imagesc(angle(Img_ref_comb(:,:,1))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(3,5,3);imagesc(angle(Img_EPI_comb(:,:,1))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(3,5,4);imagesc(angle(I_PhaseEst_LPF(:,:,1,1))); title('Phase Est after LPF'); %axis equal;square; %colormap gray;
    subplot(3,5,5);imagesc(angle(I_training_slice_hybridPh_comb(:,:,1))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
    % Slice 2 
%     subplot(3,5,6);imagesc(angle(I_PhaseEst(:,:,1,2))); title('Phase Est'); %axis equal;square; %colormap gray;
    subplot(3,5,7);imagesc(angle(Img_ref_comb(:,:,2))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(3,5,8);imagesc(angle(Img_EPI_comb(:,:,2))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(3,5,9);imagesc(angle(I_PhaseEst_LPF(:,:,1,2))); title('Phase Est after LPF'); %axis equal;square; %colormap gray;
    subplot(3,5,10);imagesc(angle(I_training_slice_hybridPh_comb(:,:,2))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
    % Slice 3
%     subplot(3,5,11);imagesc(angle(I_PhaseEst(:,:,1,3))); title('Phase Est'); %axis equal;square; %colormap gray;
    subplot(3,5,12);imagesc(angle(Img_ref_comb(:,:,3))); title('Ref comb'); %axis equal;square; %colormap gray;
    subplot(3,5,13);imagesc(angle(Img_EPI_comb(:,:,3))); title('EPI comb'); %axis equal;square; %colormap gray;
    subplot(3,5,14);imagesc(angle(I_PhaseEst_LPF(:,:,1,3))); title('Phase Est after LPF'); %axis equal;square; %colormap gray;
    subplot(3,5,15);imagesc(angle(I_training_slice_hybridPh_comb(:,:,3))); title('Ref Phase after LPF Est'); %axis equal;square; %colormap gray;
end 
% file_name = ['./Dir64E/BPU_A_GRAPPA_b',num2str(iReps),'_Slice',num2str(StartSlice),'.mat'];
% save(file_name,'I_PhaseEst','I_SMS_comb','I_GRAPPA_comb','I_PhaseEst_LPF','I_training_slice_hybridPh_comb');
end
end