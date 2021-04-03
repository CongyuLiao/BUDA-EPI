% Step3: SMS-BUDA reconstruction for distortion-free EPI
% authors: Congyu Liao, PhD, cyliao@stanford.edu
%          Xiaozhi Cao, PhD, xiaozhic@stanford.edu
% 10/20/2019



clear;clc;


addpath tools
addpath(genpath('./tools/matlab_freesurfer'))
addpath(genpath('./tools/Recon_Jon/mrir_toolbox'))
addpath svdandpca
addpath mtimesx_20110223


% parallel computing flag:
show_mercy = 10; %how many cores you wanna save for your dear colleagues, show mercy

% save path
save_path = './recon/';


% load protocol and sens map
load([save_path,'prot.mat']);
nRF = 1;  % --1 conventional EPI
          %--5 gSlider factor
nDif =34;  % 68 diffusion directions
load ([save_path,'sens_gre.mat']);
[N(1), N(2), num_chan, num_slc]= size(sens_gre);


esp = prot.iEffectiveEpiEchoSpacing * 1e-6;   % echo spacing;
AccZ = 3;                          % sms factor
AccY = prot.lAccelFactPE;     % PAT factor

% BUDA settings
opt.step_size = 0.8;
opt.num_iter = 150;   % 100 iterations
opt.tol = 0.25;
opt.winSize = [7,7]; % kernel size
opt.lambda_msl = 1.00;  % lambda of henkel low rank matrix
opt.esp=esp;   % echo-spacing
opt.PhaseShiftBase = pi;  % CAIPI shift of CAIPI
opt.show_mercy = show_mercy ;  % show how many cores you wanna use
opt.fista=1 ;  % --1 use FISTA algorithm  --0 use POCS algorithm
opt.t_k=1;



for ii_dif=1:nDif
    for ii_rf=1:nRF
        % load raw kspace data
        fpkData = [save_path,'k_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.mat'];
        load(fpkData);
        k_ap = zeros(N(1),N(2),num_chan,size(kspace_cor_ap,4));
        k_ap(:,ky_idx_ap,:,:)= kspace_cor_ap;
        k_pa = zeros(N(1),N(2),num_chan,size(kspace_cor_pa,4));
        k_pa(:,ky_idx_pa,:,:)= kspace_cor_pa;
        
        % load B0 maps
        fpField = [save_path, 'fieldmap_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii.gz'];
        dt = load_nifti(fpField);
        img_fieldmap = dt.vol;
        img_fieldmap =crop (img_fieldmap, [N(1),N(2),num_slc ]);
        
        
        tic
        % BUDA recon w/o parallel computing --it's slow, you can see the iterations
%         [img_msl]=recon_SMS_BUDA_XCLv2_fista (k_ap, k_pa,ky_idx_ap,ky_idx_pa,sens_gre, AccY, AccZ,img_fieldmap,opt);
        
        % BUDA_recon with parallel computing 
        [img_msl]=recon_SMS_BUDA_XCLv2_parfor_fista (k_ap, k_pa,ky_idx_ap,ky_idx_pa,sens_gre, AccY, AccZ,img_fieldmap,opt);
        toc
        
        apodization_para = 0.2;
        if apodization_para > 0
            kcomb = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img_msl));
            kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;
            img_msl = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
        end
        img_msl = single(img_msl);
        save([save_path, 'img_buda_rf', num2str(ii_rf),'_dif', num2str(ii_dif),'.mat'], 'img_msl')
        clear img_msl
    end
end