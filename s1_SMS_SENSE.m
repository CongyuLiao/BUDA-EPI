% Step1: SMS-SENSE reconstruction for ap and pa shots
% authors: Congyu Liao, PhD, cyliao@stanford.edu
%          Xiaozhi Cao, PhD, xiaozhic@stanford.edu
% 10/20/2019



clear;clc;close all;


addpath (genpath('./tools/'))
addpath (genpath('./Nifti_Analyze'))
addpath svdandpca
addpath mtimesx_20110223
addpath(genpath('./tools/read_meas_dat__20140924112147/'))



% parallel computing flag:
show_mercy = 25; %how many cores you wanna save for your dear colleagues, show mercy
Ncc= 16;  % number of coil compression : from 32->16 ch


% save path
save_path = './recon/';
mkdir(save_path);


%data path
file_path = './data/';
file_name = 'meas_MID00913_FID12373_BUDA_1p25inplane_R2MB3s2_ky10_b1k_34dir.dat';  % raw data of BUDA-EPI
evp_name = 'meas_MID01329_FID12797_evp_baseline_for_PAT2MB3_with_ICE_on.dat';  % scan protocols
gre_name = 'meas_MID00909_FID12369_gre_sens_2mm57slc.dat';  % gre scan for sensitivity map


f_DiffData = [file_path, file_name];
f_EvpData =[file_path, evp_name];

[~, meas.evp] = read_meas_prot(f_EvpData);  % load evp file from "evp data"
[meas.prot, ~] = read_meas_prot(f_DiffData); % load protocol from "diffusion data"


kyshift = [1,0];  % kyshift: AP =1 PA=0;
nRF = 1;  %--1 conventioanl EPI
          %--5 gSlider encodings 

AccZ = 3;                          % sms factor
AccY = meas.prot.lAccelFactPE;     % in-plane PAT factor


num_bvalue = length(meas.prot.sDiffusion_alBValue);
nDif = 34 % number of diffusion directions
meas.evp.NRepMeas =nDif*2 %(meas.prot.sDiffusion_lDiffDirections*(num_bvalue-1)+1)*(nRF*2); % number of TRs with ap and pa shots
meas.evp.RawRep =nDif*2 %(meas.prot.sDiffusion_lDiffDirections*(num_bvalue-1)+1)*(nRF*2);% number of TRs with ap and pa shots
meas.evp.RawSlc = meas.prot.sSliceArray_lSize;
meas.evp.NSlcMeas = meas.prot.sSliceArray_lSize;
evp=meas.evp;
prot =meas.prot;


mkdir(save_path);
save([save_path,'evp.mat'],'evp');
save([save_path,'prot.mat'],'prot');

PhaseShiftBase = pi; % blipped-CAIPI shift: 0- no capi shfit; pi-> FOV/2 shift 

voxel_size = [prot.dReadoutFOV, prot.dPhaseFOV, prot.dThickness] ./ ...
    [prot.lBaseResolution, prot.lPhaseEncodingLines, prot.sSliceArray_lSize]
esp = meas.prot.iEffectiveEpiEchoSpacing * 1e-6   % echo spacing
start_line = meas.evp.NFirstLin;  % the start line of EPI acquisition
N = [prot.lBaseResolution, prot.lPhaseEncodingLines]  % image size

%% load GRE data and calculate sensitivity maps
if (~exist (strcat(save_path,'sens_gre.mat'),'file'))
    fGRE=[file_path, gre_name];
    [Img_gre, whmcc]= load_GRE_data_prewht(fGRE,N,Ncc);  % load ref data, do coil compression and noise prewhitening
    % estimate ESPIRiT coil senstivity 
    sens_gre = sq(CoilSense_ESPIRIT3d( permute(Img_gre, [1,2,4,3]))); % calculate coil sensitivity using ESPIRIT
    
    sens_gre= permute(sens_gre, [1,2,4,3]);
    sens_gre= single(squeeze(sens_gre));
    save (strcat(save_path,'sens_gre.mat'), 'sens_gre','whmcc', '-V7.3');
elseif (~exist (strcat(save_path,'sens_gre'),'var'))
    load([save_path, 'sens_gre'])
end
        %%

for ii_dif=1:nDif      % diffusion directions 1:68
    for ii_rf=1:nRF    % --1 conventional EPI --5 gSlider-EPI
        for ii_rep=1:2  % two TRs:AP/PA  shots
            iReps = ii_rep+(ii_rf-1)*2+(ii_dif-1)*(2*nRF)  % odd number--> ap, even number--> pa
            % load one TR of the raw data
            [meas] = load_BUDA_SMS_data(f_DiffData,iReps,AccZ,evp);
            num_slc = size(meas.data,10)*AccZ;
            
            % ghost correction for epi data
            [~,kspace_cor] = ghost_correct_v2_STD(prot, meas.data, meas.data_phascor1d, AccY,start_line);
            
            % noise prewhitening +  coil compression
            [kspace_cor] = coilcomp_prewht(kspace_cor,whmcc);
            
            % caipirinha_deblur
            [kspace_cor]= caipi_deblur_CLv1(kspace_cor, prot, evp, AccZ,PhaseShiftBase);
            
            % flip AP and PA shots and zero-pad the partial Fourier part
            [kspace_cor]= BUDA_preprocess(kspace_cor,prot,ii_rep,kyshift);
            
            % detect ky lines of each shot
            if ii_rep ==1
                [kspace_cor_ap,ky_idx_ap]=detect_kyLines(kspace_cor); % ap shot
            else
                [kspace_cor_pa,ky_idx_pa]=detect_kyLines(kspace_cor); % pa shot
            end
        end
        [N(1), N(2), num_chan,~] = size(kspace_cor_ap);

        
        % do sms-sense for each shot
         delete(gcp('nocreate'))  
        tic
        [img_recon_pad] = recon_SMS_data_XCLv2_parfor(kspace_cor_ap, kspace_cor_pa,ky_idx_ap,ky_idx_pa,sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
        toc
        
        % zero-pad the data along slice dimension so that the number of slice can be divisible by 8
        % the zero-padded images are used for FSL topup
        num_slc_pad=ceil(num_slc/8)*8;
        img_recon_pad = zpad(img_recon_pad,[size(img_recon_pad,1),size(img_recon_pad,2),num_slc_pad,size(img_recon_pad,4)]);
        
        
        fpB0Topup = [save_path, 'img_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'];
        genNii(single( abs( (img_recon_pad)) ), voxel_size, fpB0Topup)
        
        
        kspace_cor_ap = kspace_cor_ap(:,ky_idx_ap,:,:); kspace_cor_pa = kspace_cor_pa(:,ky_idx_pa,:,:);
        kspace_cor_ap =single(kspace_cor_ap); kspace_cor_pa =single(kspace_cor_pa);
        save (strcat(save_path,'k_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.mat'), 'kspace_cor_ap', 'kspace_cor_pa','ky_idx_ap','ky_idx_pa', '-V7.3');
        
    end 
end
