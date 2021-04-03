%% ________________________________________________________________________
%% Slice and inplane GRAPPA with virtual conjugate coil reconstruction
%% Congyu Liao, PhD,  cliao2@mgh.harvard.edu
%% Martinos Center, Massachusetts General Hospital, Harvard Medical School
%% $10/08/2018
%% ______________________________________________________________________
clear all;
close all;
clc;

addpath utils


%% 
for iReps = 1:1
    AccZ = 2
    AccY = 4;
    removeROOS=0;
    extraACS=0 ;
    doMussels=0;
    do_noise_prew=1;
    if (iReps<6 || iReps >325)
        FilterSize = [40,40];
    else
        FilterSize = [30,30];
    end
    folderPath='./DATA/';
    AcqTime='20181230_nprew/';
    mkdir(strcat(folderPath,AcqTime));
    savePath1=strcat(folderPath,AcqTime,'/'); 

%     
%        PAT4 with dynamic shimming
    f_DiffData='/autofs/cluster/kawin/Congyu/CODE/SMS_Diff/raw_data/20181230/meas_MID00669_FID76973_gslider_sms2pat4_p92iso_AP_shimon.dat'
    
    if (extraACS==1)
%         f_ACSdata='/autofs/space/lucifer_002/users/congyu/20181012_shimGslider/meas_MID01346_FID55281_gSlider_1mm_PAT4SMS2_sAP_175degree_matchTETR.dat';
    end
    
    PAT_acqMethod='std';%std % FLEET
 
    Nreps = 1;
    BeginRep = iReps;
    if AccZ==1
        SMSdata=0;
    else
        SMSdata = 1;
    end
    DataIsVD = 1;
    

    [meas] = read_meas_dat_memmap_v3_1(f_DiffData,Nreps,BeginRep,SMSdata,0,DataIsVD,AccZ);
    
    if (extraACS==1)
       [meas_ACS] = read_meas_dat_memmap_v3_1(f_ACSdata,Nreps,1,SMSdata,0,DataIsVD,AccZ);   
    end
    
    if (removeROOS==0)
        if (mod (AccZ,2)==0)
            EPI_data = cat(10, meas.data(:,:,:,:,:,:,:,:,:,end),meas.data(:,:,:,:,:,:,:,:,:,1:end-1));
            % Nav of Data
            EPI_nav = cat(10,meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end),meas.data_phascor1d(:,:,:,:,:,:,:,:,:,1:end-1));
        else   
            
            EPI_data = meas.data(:,:,:,:,:,:,:,:,:,:);
            % Nav of Data
            EPI_nav = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,:);
        end         
            % Pat Ref
            if (extraACS==1)
               PAT_ref = meas_ACS.patrefscan(:,:,:,:,:,:,:,:,:,:);
            % Nav of PAT Ref
               PAT_nav = meas_ACS.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:); 
            else
               PAT_ref = meas.patrefscan(:,:,:,:,:,:,:,:,:,:);
            % Nav of PAT Ref
               PAT_nav = meas.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:); 
            end
            
            if AccZ>1
                if (extraACS==1)
                    % SMS Ref
                    SMS_ref = meas_ACS.smsrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of SMS Ref
                    SMS_nav = meas_ACS.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                else
                    % SMS Ref
                    SMS_ref = meas.smsrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of SMS Ref
                    SMS_nav = meas.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                end
            end
             
    else 
        if (mod (AccZ,2)==0)
            EPI_data = cat(10, meas.data(:,:,:,:,:,:,:,:,:,end),meas.data(:,:,:,:,:,:,:,:,:,1:end-1));
            % Nav of Data
            EPI_nav = cat(10,meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end),meas.data_phascor1d(:,:,:,:,:,:,:,:,:,1:end-1));
        else   
            
            EPI_data = meas.data(:,:,:,:,:,:,:,:,:,:);
            % Nav of Data
            EPI_nav = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,:);
        end
           % Nav of Data
            EPI_data =removeReadoutOS (EPI_data);
            EPI_nav = removeReadoutOS (EPI_nav);
            
                if (extraACS==1)
                    % Pat Ref
                    PAT_ref = meas_ACS.patrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of PAT Ref
                    PAT_nav = meas_ACS.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                    
                    PAT_ref = removeReadoutOS (PAT_ref);
                    PAT_nav = removeReadoutOS (PAT_nav);
                else                    
                    % Pat Ref
                    PAT_ref = meas.patrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of PAT Ref
                    PAT_nav = meas.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                    
                    PAT_ref = removeReadoutOS (PAT_ref);
                    PAT_nav = removeReadoutOS (PAT_nav);
                end
           
            if AccZ>1
                if (extraACS==1)
                    % SMS Ref
                    SMS_ref = meas_ACS.smsrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of SMS Ref
                    SMS_nav = meas_ACS.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                    
                    SMS_ref = removeReadoutOS (SMS_ref);
                    SMS_nav = removeReadoutOS (SMS_nav);
                else
                    % SMS Ref
                    SMS_ref = meas.smsrefscan(:,:,:,:,:,:,:,:,:,:);
                    % Nav of SMS Ref
                    SMS_nav = meas.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
                    
                    SMS_ref = removeReadoutOS (SMS_ref);
                    SMS_nav = removeReadoutOS (SMS_nav);
                end
            end
 
    end
    kSize_GRAPPA = [5,5];  % grappa kernel
    lambda_tik= eps;  % lambda for grappa regularization
    num_acs = [200,100];          % ACS size 

    nSlice = size(meas.data,10)*AccZ;
    
    if (removeROOS==0)
        nLines = size(meas.data,1);
    else
        nLines = size(meas.data,1)/2;
    end
    nColumnData = size(meas.data,2);
    nColumn = size(meas.data,2)*AccY;
    PhaseShiftBase = pi;  %PE shift 1-pi 2-pi/2
        % to find startline of EPI data
    if AccZ>1
        for sL=1:AccY
            tmp=squeeze(SMS_ref);
            tmp1=tmp(:,sL,:,:,:);
            tmp2=sum(tmp1(:));
            if tmp2~=0
                StartLine = sL
            end
            clear tmp tmp1 tmp2
        end
    else
        StartLine = 2;
    end
        
    clear sL
    %coil compression
    tailor_SMS_ghost = 0;
    process_REF_Nav = 1;
    process_PAT_GHOST_LOAD = 1;

    nLines_old = nLines;
    NRead = [(nLines/2-95+1):(nLines/2+95)];
    nLines = length(NRead);
    NormCoils = [1:20];  % 64     
    nCoil = length(NormCoils);
       
    if nCoil== size(meas.data,3);
        coil_comp =0  %0-- without coil compression
    else
        coil_comp =1  %1-- with coil compression  
    end
    if coil_comp == 1
       if iReps == 1
            img_to_compress = permute(squeeze(sum(EPI_data(:,:,:,:,1,:,1,:,1,:),8)),[1 2 4 3]);
            [ ~, comp_mtx]  =  svd_compress3d( img_to_compress, nCoil, 0 );
            save(strcat(savePath1,'comp_mtx.mat'), 'comp_mtx', '-V7.3')
       end
       load(strcat(savePath1,'comp_mtx.mat'))
       EPI_data= svd_apply_all(EPI_data,comp_mtx);
       EPI_nav = svd_apply_all(EPI_nav ,comp_mtx);
       PAT_ref = svd_apply_all(PAT_ref ,comp_mtx);
       PAT_nav = svd_apply_all(PAT_nav ,comp_mtx);
       if AccZ>1
           SMS_ref = svd_apply_all(SMS_ref ,comp_mtx);
           SMS_nav = svd_apply_all(SMS_nav ,comp_mtx);
       end
    
    end
    % phase correction for PAT references
    if ( exist (strcat(savePath1,'kspace_acs_all.mat'),'file'))    
        % Load Corrected PAT Ref from a file
        load (strcat(savePath1,'kspace_acs_all.mat') );
    else
        if (strcmp(PAT_acqMethod,'FLEET')==1)
            PAT_ref_PATnavCor_img  = ghost_correct_pat_ref_v1_FLEET(meas.prot, PAT_ref,PAT_nav, AccY);
        
        else
           PAT_ref_PATnavCor_img  = ghost_correct_pat_ref_v1_STD(meas.prot, PAT_ref,PAT_nav, AccY);
        end
        
        % reduce numbers of coils
        temp = PAT_ref_PATnavCor_img;PAT_ref_PATnavCor_img = [];
        PAT_ref_PATnavCor_img = temp(:,:,NormCoils,1,1,1,1,:,1,:);
        clear temp;

        % Prepare training data for in-plane GRAPPA
        Krecon_PAT_ref_PATnavCor_mslc = mrir_fDFT_freqencode(mrir_fDFT_phasencode(PAT_ref_PATnavCor_img(:,:,:,1,1,1,1,1,1,:)));
        Krecon_PAT_ref_PATnavCor_mslc_ext = zpad(squeeze(Krecon_PAT_ref_PATnavCor_mslc),nLines_old,nColumn,nCoil,nSlice);
        Krecon_PAT_ref_PATnavCor_mslc_acs_ext = zeros(nLines_old,nColumn,nCoil,1,1,1,1,1,1,nSlice);
        Krecon_PAT_ref_PATnavCor_mslc_acs_ext(:,:,:,1,1,1,1,1,1,:) = Krecon_PAT_ref_PATnavCor_mslc_ext;
        Irecon_PAT_ref_PATnavCor_mslc_acs_ext = mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_PAT_ref_PATnavCor_mslc_acs_ext));
        img_PAT_acs = squeeze(Irecon_PAT_ref_PATnavCor_mslc_acs_ext (:,:,:,:));
        k_PAT_acs=squeeze(Krecon_PAT_ref_PATnavCor_mslc_acs_ext);
        kspace_acs_all=zeros(size(k_PAT_acs));
        kspace_acs_all(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2 + 1, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2 + 1,:, :) = k_PAT_acs(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2 + 1, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2 + 1, :,:);

        clear Krecon_PAT_ref_PATnavCor_mslc  Krecon_PAT_ref_PATnavCor_mslc_ext Krecon_PAT_ref_PATnavCor_mslc_acs_ext
        clear Irecon_PAT_ref_PATnavCor_mslc_acs_ext
        save(strcat(savePath1,'kspace_acs_all.mat'), 'kspace_acs_all', '-V7.3')
   
    end 
    
    
    if (do_noise_prew)  
         if ( exist (strcat(savePath1,'noise_prew.mat'),'file'))  
             load (strcat(savePath1,'noise_prew.mat') );
         else          
            ref=img_PAT_acs(:,:,:,floor(nSlice/2));
            [psi,L_inv,pos]=estimate_psi(ref,10);
            save(strcat(savePath1,'noise_prew.mat'), 'psi', 'L_inv','pos','-V7.3')
         end
        kspace_acs_all=data_nprew(kspace_acs_all,L_inv);
        img_PAT_acs=data_nprew(mrir_iDFT_freqencode(mrir_iDFT_phasencode(kspace_acs_all)),L_inv);
        EPI_data= data_nprew2(EPI_data,L_inv);
        EPI_nav = data_nprew2(EPI_nav ,L_inv);
        PAT_ref = data_nprew2(PAT_ref ,L_inv);
        PAT_nav = data_nprew2(PAT_nav ,L_inv);
        if AccZ>1
           SMS_ref = data_nprew2(SMS_ref ,L_inv);
           SMS_nav = data_nprew2(SMS_nav ,L_inv);
        end
    end
    
    
     % sensitivity maps
    if (exist (strcat(savePath1,'SensMap_PAT_ref.mat'),'file'))
        load(strcat(savePath1,'SensMap_PAT_ref.mat'))
    else
         SensMap_all=CoilSense_ESPIRIT3d(permute(img_PAT_acs,[1 2 4 3]));
         SensMap_all=permute(SensMap_all,[1 2 4 3]);
         save (strcat(savePath1,'SensMap_PAT_ref.mat'), 'SensMap_all', '-V7.3');
    end
    if AccZ>1
        for zz = 1:nSlice
            DThickness(zz) = meas.prot.sSliceArray(zz).dThickness;
        end
    end
 %% Step 1: Conventional slice- and inplane-GRAPPA reconstruction    
    matlabpool 4
parfor iSlice=1:17
% for iSlice=3:3

    AccZ = 2;
    AccY = 4;
    StartSlice =iSlice;
    SensMap= SensMap_all(:,:,:,StartSlice:(nSlice/AccZ):nSlice);
    kspace_acs= kspace_acs_all(:,:,:,StartSlice:(nSlice/AccZ):nSlice);
    file_name=strcat(PAT_acqMethod,'_gslider_EPI_SMS',num2str(AccZ),'PAT',num2str(AccY),'_sSlice',num2str(StartSlice));
    mkdir(strcat(folderPath,AcqTime),file_name);
    savePath=strcat(folderPath,AcqTime,file_name,'/'); 

    
    TYPE_RECON = 'VCC-GRAPPA';
%     TYPE_CalGFactMap = 'G-FACTOR'; 
    TYPE_CalGFactMap = 'no_gfactor calculation' ;
    TYPE_SAVE = 'SAVE';

    NReps = iReps;
    disp(['StartSlice is ',num2str(StartSlice)])
    disp(['Number of Reps is ',num2str(NReps)])
    
    %  EPI ghost corrections
    lin_fit_EPInavCor = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav);
    EPI_data_EPInavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(EPI_data, lin_fit_EPInavCor);
    data_hybrid = mrir_iDFT_freqencode(EPI_data_EPInavCor); 
    EPI_data_EPInavCor=[];
    data_hybrid_tp = mrir_regrid_trapezoid(data_hybrid, meas.prot);
    EPI_data_EPInavCor = mrir_fDFT_freqencode(data_hybrid_tp);
%     clear data_hybrid data_hybrid_tp
    data_hybrid=[];data_hybrid_tp=[];
    
 if (AccZ>1)
    dims_EPI_nav = size(EPI_nav);
    dims_SMS_nav = size(SMS_nav);
    AccZ1 = dims_SMS_nav(10)/dims_EPI_nav(10);
    if (AccZ1~=AccZ)
        disp(['slice acceleration factor is ', num2str(AccZ1)])
        AccZ=AccZ1;
    end
%     clear AccZ1
    AccZ1=[];
    dims_EPI_data = size(EPI_data);
    dims_SMS_ref = size(SMS_ref);
    if tailor_SMS_ghost == 1
        if length(dims_EPI_nav) < 10
            dims_EPI_nav(end+1:10) = 1; 
        end
        lin_fit_EPInavCor = reshape( median( reshape(lin_fit_EPInavCor, [2 1 dims_EPI_nav(3:7) 1 1 dims_EPI_nav(10)]), 7) , 2,[]) ; % for saving to use for tailor ghost correction   
    end
    lin_fit_EPInavCor_reshape = reshape(lin_fit_EPInavCor,2,dims_EPI_nav(3),dims_EPI_nav(7),dims_EPI_nav(10));
    lin_fit_EPInavCor_SMSref_reshape = mean(repmat(lin_fit_EPInavCor_reshape,[1,1,1,AccZ]),3);
    lin_fit_EPInavCor_SMSref= reshape(lin_fit_EPInavCor_SMSref_reshape,2,dims_SMS_nav(3)*dims_SMS_nav(10));
     lin_fit_EPInavCor=[]; lin_fit_EPInavCor_reshape =[];
     lin_fit_EPInavCor_SMSref_reshape=[];
    dims_EPI_nav=[]; dims_SMS_nav=[]; dims_EPI_data=[]; dims_SMS_ref=[];

    % averaged ghost correction by using SMS navigators
    SMS_ref_EPInavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(SMS_ref, lin_fit_EPInavCor_SMSref);
    smsrefscan_hybrid = mrir_iDFT_freqencode(SMS_ref_EPInavCor); 
	SMS_ref_EPInavCor = [];
    smsrefscan_hybrid_tp = mrir_regrid_trapezoid(smsrefscan_hybrid, meas.prot);
    SMS_ref_EPInavCor = mrir_fDFT_freqencode(smsrefscan_hybrid_tp);
    smsrefscan_hybrid=[]; smsrefscan_hybrid_tp=[];

    % Slice-GRAPPA reconstruction
    % CaipiShift ACS data 
    temp = sum(SMS_ref_EPInavCor(:,:,:,:,:,:,1,:,:,:),8);
    K_SG_acs = temp(:,StartLine:AccY:end,:,:,:,:,1,:,:,:);
     temp=[];
    K_SG_acs_mslc = K_SG_acs(:,:,:,:,:,:,1,:,:,StartSlice:(nSlice/AccZ):nSlice);
    SliceGroup = [0:(AccZ-1)]
    K_SG_acs_mslc_shft = zeros(size(K_SG_acs_mslc));
    for ss = 1:AccZ
        K_SG_acs_mslc_shft(:,:,:,:,:,:,1,:,:,ss) = CaipirinhaShift_K_v2(K_SG_acs_mslc(:,:,:,:,:,:,1,:,:,ss),SliceGroup(ss),PhaseShiftBase);
    end
    I_SG_acs_mslc_shft= mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_acs_mslc_shft));
    SliceGroup=[]; K_SG_acs_mslc_shft=[];  K_SG_acs_mslc=[]; K_SG_acs=[];
    
    %  CaipiShift EPI data 
    K_SG_EPI = sum(EPI_data_EPInavCor(:,:,:,:,:,:,1,:,:,:),8); % 1st TR
    
    SliceSep = sum(DThickness(StartSlice:(StartSlice+(nSlice/AccZ)-1)));
    %SliceSep = (nSlice/AccZ)*DThickness(1);
    K_SG_EPI_deblur = CaipirinhaDeblur_v4(K_SG_EPI, meas.prot, meas.evp, PhaseShiftBase, SliceSep);
    K_SG_EPI_1slc_deblur =  K_SG_EPI_deblur(:,:,:,:,:,:,1,:,:,StartSlice);
    I_SG_EPI_1slc_deblur = mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_EPI_1slc_deblur));
    K_SG_EPI_deblur=[]; K_SG_EPI_1slc_deblur=[]; K_SG_EPI =[];SliceSep=[];
    
    %  slice-GRAPPA recon
    disp(['Slice:', num2str(StartSlice)]);
    img_EPI= permute(squeeze(I_SG_EPI_1slc_deblur),[3 2 1]);
    img_acs = permute( squeeze(I_SG_acs_mslc_shft),[3 2 1 4]);
    [Irecon_SG,Irecon_SG_sep,gfactor,weights1_Slice,weights2_Slice] = SliceGRAPPA_v6_3_wo_lambda(img_EPI,img_acs,3,3,PhaseShiftBase);%% STD
     img_EPI=[]; img_acs=[];
    
    %  correct recon data of slice-GRAPPA with the residul of Navigates of SMS ref and Data
    % prepare the kspace data +mask even/odd lines
    temp = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    temp(:,:,:,1,1,1,1,1,1,:) = permute(Irecon_SG_sep(:,:,:,:),[3 2 1 4]);
    Krecon_SG_sep = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Krecon_SG_sep  = mrir_fDFT_freqencode(mrir_fDFT_phasencode(temp));
    temp=[];
    
    Krecon_SG_sep_mask = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,2,1,AccZ);
    Mask_odd = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_odd(:,2:2:end,:,:) = 0;
    Mask_even = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_even(:,1:2:end,:,:) = 0;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,1,1,:) =Krecon_SG_sep .*Mask_odd;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,2,1,:) =Krecon_SG_sep .*Mask_even;
    Krecon_SG_sep =[]; Mask_odd=[]; Mask_even=[];
    
    %  add EPI_nav back (inversion)
    EPI_nav1slc = EPI_nav(:,:,:,:,:,:,1,:,:,StartSlice);
    lin_fit_EPInavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav1slc);
    Krecon_SG_sep_EPInavCor_inv= zeros(nLines_old,nColumnData,nCoil,1,1,1,1,2,1,AccZ);
    for ss = 1:AccZ
        Krecon_SG_sep_EPInavCor_inv(:,:,:,1,1,1,1,:,1,ss) = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_mask(:,:,:,1,1,1,1,:,1,ss), -lin_fit_EPInavCor1slc );
    end
       
    % correction using SMS_nav 
	SMS_nav1slc= SMS_nav(:,:,:,1,1,1,1,:,1,StartSlice:(nSlice/AccZ):nSlice);
    lin_fit_SMSnavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav1slc);
    temp_corr = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_EPInavCor_inv, lin_fit_SMSnavCor1slc);
    Krecon_SG_sep_SMSnavCor = sum(temp_corr(:,:,:,1,1,1,1,:,1,:),8);
    Irecon_SG_sep_SMSnavCor = squeeze(mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_SG_sep_SMSnavCor)));
    temp_corr=[]; SMS_nav1slc=[]; lin_fit_SMSnavCor1slc=[]; EPI_nav1slc=[]; lin_fit_EPInavCor1slc=[]; 
    Krecon_SG_sep_EPInavCor_inv=[];
    temp = squeeze(sos(Irecon_SG_sep_SMSnavCor(:,:,:,2),3,2)-sos(permute(Irecon_SG_sep(:,:,:,2),[3 2 1]),3,2));
    
%     mosaic(squeeze(sos(Irecon_SG_sep_SMSnavCor(:,:,:,2),3,2)-sos(permute(Irecon_SG_sep(:,:,:,2),[3 2 1]),3,2)),1,1,41,['mean of residul = ',num2str(mean_delta_residul)]);
    temp=[]; mean_delta_residul=[];
 
    %get ghost-free SMS_ref data
    
    lin_fit_SMSnavCor= mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav);
	SMS_ref_SMSnavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(SMS_ref,lin_fit_SMSnavCor);
    smsrefscan_hybrid = mrir_iDFT_freqencode(SMS_ref_SMSnavCor); 
    SMS_ref_SMSnavCor = [];
    smsrefscan_free_hybrid = mrir_regrid_trapezoid(smsrefscan_hybrid, meas.prot);
    SMS_ref_SMSnavCor= mrir_fDFT_freqencode(smsrefscan_free_hybrid);
     smsrefscan_hybrid=[]; smsrefscan_free_hybrid =[];lin_fit_SMSnavCor=[];
  
 
    temp = sum(SMS_ref_SMSnavCor(:,:,:,:,:,:,1,:,:,:),8);
    Krecon_SMS_ref_SMSnavCor = temp(:,2:AccY:end,:,:,:,:,1,:,:,:);
    temp=[];
    Krecon_sms_ref_mslc = Krecon_SMS_ref_SMSnavCor(:,:,:,:,:,:,1,:,:,StartSlice:(nSlice/AccZ):nSlice);
    Irecon_sms_ref_mslc = mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_sms_ref_mslc));
    Irecon_sms_ref_mslc = permute( squeeze(Irecon_sms_ref_mslc),[3 2 1 4]);
 end
    % in-plane GRAPPA reconstruction
    %  in-plane GRAPPA Recon
    Irecon_EPI_Gfree_fs = zeros(nLines_old,nColumn,nCoil,AccZ);
    Irecon_PAT_acs = zeros(nLines_old,nColumn,nCoil,AccZ);
    
    Irecon_EPI_Gfree_fs_comb = zeros(nLines_old,nColumn,AccZ);
    Irecon_PAT_acs_comb = zeros(nLines_old,nColumn,AccZ);
    
    ksampled_EPI=zeros(nLines_old,nColumn,nCoil,AccZ);
    if AccZ>1
        Irecon_SMS_ref_fs = zeros(nLines_old,nColumn,nCoil,AccZ);
        Irecon_SMS_ref_fs_comb = zeros(nLines_old,nColumn,AccZ);
        ksampled_EPI(:,StartLine:AccY:end,:,:)=squeeze(Krecon_SG_sep_SMSnavCor);
    else
        tmp=squeeze(sum(EPI_data_EPInavCor(:,:,:,:,:,:,1,:,:,:),8));
        ksampled_EPI(:,StartLine:AccY:end,:,:)=squeeze(tmp(:,:,:,StartSlice));
        tmp=[];
    end
    if AccZ>1
        ksampled_SMS=zeros(nLines_old,nColumn,nCoil,AccZ);
        ksampled_SMS(:,StartLine:AccY:end,:,:)=squeeze(Krecon_sms_ref_mslc);
    end
    weights_Inplane=zeros(nCoil, nCoil,nColumn,nLines_old,AccZ );
    for slc = 1:AccZ
    % in-plane GRAPPA for ky undersampled, ghost-free EPI data (Irecon_SG_sep_SMSnavCor) 
   
    [Irecon_EPI_Gfree_fs(:,:,:,slc),img_weights] = grappa_gfactor_2d_vc_CL_v2( ksampled_EPI(:,:,:,slc),kspace_acs(:,:,:,slc), 1, AccY, num_acs, kSize_GRAPPA, lambda_tik ,0, zeros(1,nCoil), zeros(1,nCoil),StartLine);  
    Irecon_EPI_Gfree_fs_comb(:,:,slc)= sum(conj(SensMap(:,:,:,slc)).*(Irecon_EPI_Gfree_fs(:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
    if AccZ>1
        [Irecon_SMS_ref_fs(:,:,:,slc),~] =grappa_gfactor_2d_vc_CL_v3(ksampled_SMS(:,:,:,slc),kspace_acs(:,:,:,slc), 1, AccY, num_acs, kSize_GRAPPA, lambda_tik,0, zeros(1,nCoil), zeros(1,nCoil) ,StartLine);
        Irecon_SMS_ref_fs_comb(:,:,slc) = sum(conj(SensMap(:,:,:,slc)).*(Irecon_SMS_ref_fs(:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
    end
    Irecon_PAT_acs (:,:,:,slc) = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kspace_acs(:,:,:,slc))) ;
    Irecon_PAT_acs_comb(:,:,slc) = sum(conj(SensMap(:,:,:,slc)).*(Irecon_PAT_acs (:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
     
    end
%     mosaic(Irecon_SMS_ref_fs_comb(NRead,:,:),1,3,6,'Recon of SMS Ref');%caxis([0 100]);
%     mosaic(Irecon_PAT_acs_comb(NRead,:,:),1,3,7,'Recon of PAT Ref');%caxis([0 100]);
%     mosaic(Irecon_EPI_Gfree_fs_comb(NRead,:,:),1,3,8,['Recon of Mag(dwi) data b',num2str(iReps)]);%caxis([0 50]);
%   

%     mosaic(Irecon_EPI_Gfree_fs_comb(:,:,:),1,3,8,['Recon of dwi data b',num2str(1)]);%caxis([0 50]);

%% Step 2: background phase estimation
    %2.1 Estimation background phase of SMS Ref------------------------------------------
 if (AccZ>1)
    disp('Step 2.1:  Estimation backbround phase from slice VCC-GRAPPA reconstruction data')
    [Irecon_SMS_bkgEst,Krecon_SMS_bkgEst,bkgPhase,PhaseDiff_LPF]=estimate_background_phase_v2 (Irecon_EPI_Gfree_fs_comb,Irecon_SMS_ref_fs_comb,Irecon_SMS_ref_fs,SensMap,FilterSize);

    %2.2 Estimation backbround phase of PAT Ref------------------------------------------
    disp('Step 2.2:  Estimation backbround phase from inplane VCC-GRAPPA reconstruction data')
    if (strcmp(PAT_acqMethod,'FLEET')==1)
        % for feet acquisition, use estimated phase from sms-ref 
        [Irecon_PAT_bkgEst,Krecon_PAT_bkgEst]=estimate_background_phase_fleet (Irecon_EPI_Gfree_fs_comb,Irecon_PAT_acs_comb,Irecon_PAT_acs,SensMap,PhaseDiff_LPF,Irecon_SMS_bkgEst);  
    else
        % for std acquisition, use its own phase
        [Irecon_PAT_bkgEst,Krecon_PAT_bkgEst,~,~]=estimate_background_phase_v2 (Irecon_EPI_Gfree_fs_comb,Irecon_PAT_acs_comb,Irecon_PAT_acs,SensMap,FilterSize);
    end
 else
     disp('Step 2.2:  Estimation backbround phase from inplane VCC-GRAPPA reconstruction data')
     [Irecon_PAT_bkgEst,Krecon_PAT_bkgEst,~,~]=estimate_background_phase_v2 (Irecon_EPI_Gfree_fs_comb,Irecon_PAT_acs_comb,Irecon_PAT_acs,SensMap,FilterSize);
 end  
%% Step 3: Virtual Coil Concept GRAPPA (VCC-GRAPPA) recon
if (AccZ>1)    
%  Slice VCC-GRAPPA
    disp('Step 3.1: VCC-Slice-GRAPPA reconstruction')
    % Prepare slice_acs

    temp = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    temp(:,:,:,1,1,1,1,1,1,:) = AccY.*Krecon_SMS_bkgEst(:,StartLine:AccY:end,:,:);
	Krecon_SMS_bkgEst_mask = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,2,1,AccZ);
    Mask_odd = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_odd(:,1:2:end,:,:) = 0;
    Mask_even = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_even(:,2:2:end,:,:) = 0;
    Krecon_SMS_bkgEst_mask(:,:,:,1,1,1,1,1,1,:) = temp.*Mask_odd;
    Krecon_SMS_bkgEst_mask(:,:,:,1,1,1,1,2,1,:) = temp.*Mask_even;
    temp =[];Mask_odd=[]; Mask_even=[]; Krecon_SMS_bkgEst=[];
    Krecon_SMS_bkgEst=Krecon_SMS_bkgEst_mask;
    Krecon_SMS_bkgEst_mask=[];
    
    %  add Nav of SMS Ref back (inversion)
    SMS_nav1slc=[];  lin_fit_SMSnavCor1slc=[];
    SMS_nav1slc = SMS_nav(:,:,:,1,1,1,1,:,1,StartSlice:(nSlice/AccZ):end);
    lin_fit_SMSnavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav1slc);
    Krecon_SMS_bkgEst_SMSnavCor_inv = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SMS_bkgEst, -lin_fit_SMSnavCor1slc);
    SMS_nav1slc=[];  lin_fit_SMSnavCor1slc=[];
    
    %correct epi data using Nav of EPI
    lin_fit_EPInavCor=[];
    lin_fit_EPInavCor = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav);
    
    dims_EPI_nav = size(EPI_nav);
    dims_SMS_nav = size(SMS_nav);   
    

    lin_fit_EPInavCor_reshape = reshape(lin_fit_EPInavCor,2,dims_EPI_nav(3),dims_EPI_nav(7),dims_EPI_nav(10));
    lin_fit_EPInavCor_SMSref_reshape = mean(repmat(lin_fit_EPInavCor_reshape,[1,1,1,AccZ]),3);
    lin_fit_EPInavCor_SMSref= reshape(lin_fit_EPInavCor_SMSref_reshape,2,dims_SMS_nav(3)*dims_SMS_nav(10));
    lin_fit_EPInavCor=[]; lin_fit_EPInavCor_reshape =[];lin_fit_EPInavCor_SMSref_reshape=[];
    dims_EPI_nav=[]; dims_SMS_nav =[];dims_EPI_data =[];dims_SMS_ref=[];
    temp = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SMS_bkgEst_SMSnavCor_inv, lin_fit_EPInavCor_SMSref);
    K_SG_acs_mslc = sum(temp,8);
    Krecon_SMS_bkgEst_SMSnavCor_inv=[]; temp=[];
    
    %CaipiShift ACS data
     SliceGroup=[]; K_SG_acs_mslc_shft=[];
    SliceGroup = [0:(AccZ-1)]
    K_SG_acs_mslc_shft = zeros(size(K_SG_acs_mslc));

    for ss = 1:AccZ
        K_SG_acs_mslc_shft(:,:,:,:,:,:,1,:,:,ss) = CaipirinhaShift_K_v2(K_SG_acs_mslc(:,:,:,:,:,:,1,:,:,ss),SliceGroup(ss),PhaseShiftBase);
    end
    I_SG_acs_mslc_shft= mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_acs_mslc_shft));
    SliceGroup=[]; K_SG_acs_mslc_shft=[];  K_SG_acs_mslc=[]; K_SG_acs=[];

   
    % CaipiShift EPI data 
    K_SG_EPI = sum(EPI_data_EPInavCor(:,:,:,:,:,:,1,:,:,:),8); % 1st TR
    
    SliceSep = sum(DThickness(StartSlice:(StartSlice+(nSlice/AccZ)-1)));
    %SliceSep = (nSlice/AccZ)*DThickness(1);
    K_SG_EPI_deblur = CaipirinhaDeblur_v4(K_SG_EPI, meas.prot, meas.evp, PhaseShiftBase, SliceSep);
    K_SG_EPI_1slc_deblur =  K_SG_EPI_deblur(:,:,:,:,:,:,1,:,:,StartSlice);
    I_SG_EPI_1slc_deblur = mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_EPI_1slc_deblur));
    K_SG_EPI_deblur=[]; K_SG_EPI_1slc_deblur=[]; K_SG_EPI=[]; SliceSep=[];
      
    %Extent VCC
    img_EPI= permute( squeeze(I_SG_EPI_1slc_deblur),[3 2 1]); 
    img_acs = permute(squeeze(I_SG_acs_mslc_shft),[3 2 1 4]);
    img_EPI_vcc = cat(1,img_EPI,conj(img_EPI));
    img_acs_vcc = cat(1,img_acs,conj(img_acs));
    img_EP=[]; Img_acs=[];

    % slice-GRAPPA recon
    [Irecon_SG_vcc,Irecon_SG_sep_vcc,gfactor,weights1_Slice_vcc,weights2_Slice_vcc] = SliceGRAPPA_v6_3_wo_lambda(img_EPI_vcc,img_acs_vcc,3,3,PhaseShiftBase);%% STD   
    img_EPI_vcc=[]; img_acs_vcc=[];



    % correct recon data of slice-VCC-GRAPPA with the residul of Nav of SMS Ref and Nav of Data
    
    % prepare the kspace data 
    Krecon_SG_sep=[];
    temp = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    temp(:,:,:,1,1,1,1,1,1,:) = permute(Irecon_SG_sep_vcc(1:nCoil,:,:,:),[3 2 1 4]);
    Krecon_SG_sep = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Krecon_SG_sep  = mrir_fDFT_freqencode(mrir_fDFT_phasencode(temp));
    temp=[];
    
    Krecon_SG_sep_mask = zeros(nLines_old,nColumnData,nCoil,1,1,1,1,2,1,AccZ);
    Mask_odd = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_odd(:,2:2:end,:,:) = 0;
    Mask_even = ones(nLines_old,nColumnData,nCoil,1,1,1,1,1,1,AccZ);
    Mask_even(:,1:2:end,:,:) = 0;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,1,1,:) =Krecon_SG_sep .*Mask_odd;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,2,1,:) =Krecon_SG_sep .*Mask_even;
    Krecon_SG_sep=[];  Mask_odd=[]; Mask_even=[];
    
    %add EPI_nav back (inversion)
    EPI_nav1slc = EPI_nav(:,:,:,:,:,:,1,:,:,StartSlice);
    lin_fit_EPInavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav1slc);
    Krecon_SG_sep_EPInavCor_inv= zeros(nLines_old,nColumnData,nCoil,1,1,1,1,2,1,AccZ);
    for ss = 1:AccZ
        Krecon_SG_sep_EPInavCor_inv(:,:,:,1,1,1,1,:,1,ss) = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_mask(:,:,:,1,1,1,1,:,1,ss), -lin_fit_EPInavCor1slc );
    end
       
    %correction using SMS_nav 
	SMS_nav1slc= SMS_nav(:,:,:,1,1,1,1,:,1,StartSlice:(nSlice/AccZ):nSlice);
    lin_fit_SMSnavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav1slc);
    temp_corr = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_EPInavCor_inv, lin_fit_SMSnavCor1slc);
    Krecon_SG_sep_SMSnavCor = sum(temp_corr(:,:,:,1,1,1,1,:,1,:),8);
    Irecon_SG_sep_SMSnavCor = squeeze(mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_SG_sep_SMSnavCor)));
    temp_corr=[]; SMS_nav1slc=[]; lin_fit_SMSnavCor1slc=[];
    EPI_nav1slc=[]; lin_fit_EPInavCor1slc =[];
     Krecon_SG_sep_EPInavCor_inv=[];
    temp = squeeze(sos(Irecon_SG_sep_SMSnavCor(:,:,:,2),3,2)-sos(permute(Irecon_SG_sep(:,:,:,2),[3 2 1]),3,2));
    mean_delta_residul = mean(temp(:));
    temp=[]; mean_delta_residul=[];
end    
  % VCC in-plane GRAPPA recon  

    % Extend VCC 
    ksampled_EPI_vcc=zeros(nLines_old,nColumn,2*nCoil,AccZ);
    ksampled_temp=zeros(nLines_old,nColumn,nCoil,AccZ);
    if (AccZ>1)     
        ksampled_temp(:,StartLine:AccY:end,:,:)=Krecon_SG_sep_SMSnavCor;
    else
        tmp=squeeze(sum(EPI_data_EPInavCor(:,:,:,:,:,:,1,:,:,:),8));
        ksampled_temp(:,StartLine:AccY:end,:,:)=squeeze(tmp(:,:,:,StartSlice));
    end
    ksampled_EPI_vcc = cat( 3, ksampled_temp, mrir_fDFT_freqencode(mrir_fDFT_phasencode(conj(mrir_iDFT_freqencode(mrir_iDFT_phasencode(ksampled_temp)) ))));
    ksampled_EPI_vcc( abs(ksampled_EPI_vcc) < 1e-12 ) = 0;
      

    % combine real and virtual coils
     temp = conj(Irecon_PAT_bkgEst); 
    img_PAT_acs_vcc = cat(3,Irecon_PAT_bkgEst,temp);%.*sqrt(2);
    k_PAT_acs_vcc = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img_PAT_acs_vcc));
    
    temp=[]; size_calib =[];DeltaY_calib=[]; DeltaY_fullSize=[];
 
    % set ksampled & kacs data
    
    kspace_acs_vcc = zeros(nLines_old,nColumn,2*nCoil,AccZ);

    
    kspace_acs_vcc(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2 + 1, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2 + 1, :,:) = ...
    k_PAT_acs_vcc(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2 + 1, 1+end/2-num_acs(2)/2:end/2+num_acs(2)/2 + 1, :,:);
    weights_Inplane_vcc=zeros(2*nCoil, 2*nCoil,nColumn,nLines_old,AccZ );
    
    Irecon_EPI_Gfree_fs_vcc=zeros(nLines_old,nColumn,2*nCoil,AccZ);
    Irecon_EPI_Gfree_fs_comb_vcc=zeros(nLines_old,nColumn,AccZ);
    % k-space offset for virtual coil: 
    DC_line=floor(nColumn/2+1);
    delY=mod(DC_line-StartLine,AccY);
    delY_vc=mod(DC_line+delY,AccY); 
    if (delY_vc==0)
        delY_vc=AccY;
    end
    shiftY_vc=delY_vc-StartLine;

    DelZ = cat(2, zeros(1, nCoil), repmat(mod(nLines_old, 1), [1, nCoil]));
    DelY = cat(2, zeros(1, nCoil), repmat(shiftY_vc, [1, nCoil]));

    %in-plane VCC-GRAPPA recon
    
    for slc = 1:AccZ
        % in-plane GRAPPA for ky undersampled, ghost-free EPI data (Irecon_SG_sep_SMSnavCor)    
        
 
        [Irecon_EPI_Gfree_fs_vcc(:,:,:,slc),img_weights_vcc] = grappa_gfactor_2d_vc_CL_v2( ksampled_EPI_vcc(:,:,:,slc),kspace_acs_vcc(:,:,:,slc), 1, AccY, num_acs, kSize_GRAPPA, lambda_tik,0, DelZ, DelY ,StartLine);      
        Irecon_EPI_Gfree_fs_comb_vcc(:,:,slc)= sum(conj(SensMap(:,:,:,slc)).*(Irecon_EPI_Gfree_fs_vcc(:,:,1:nCoil,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
    end
 %%   
%     if (doMussels==1)
%         % odd-even MUSSELS + FISTA: hard / soft thresholding  
%         for slc =1:AccZ
% %             for slc=2
%             kspace_tmp=ksampled_EPI_vcc(:,:,:,slc);
%             kspace_slice =zeros([size(kspace_tmp),2]) ;
%             kspace_slice(:,StartLine:AccY*2:end,:,1) =kspace_tmp(:,StartLine:AccY*2:end,:);
%             kspace_slice(:,StartLine+AccY:AccY*2:end,:,2) =kspace_tmp(:,StartLine+AccY:AccY*2:end,:);
%             kspace_tmp=[];
%             img_tmp=zeros(size(kspace_slice));
%       
% %             [img_tmp(:,:,:,1)] = grappa_gfactor_2d_vc_CL_v2( kspace_slice(:,:,:,1),kspace_acs_vcc(:,:,:,slc), 1, AccY*2, num_acs, kSize_GRAPPA, lambda_tik,0, DelZ, DelY ,StartLine);      
% %             [img_tmp(:,:,:,2)] = grappa_gfactor_2d_vc_CL_v2( kspace_slice(:,:,:,2),kspace_acs_vcc(:,:,:,slc), 1, AccY*2, num_acs, kSize_GRAPPA, lambda_tik,0, DelZ, DelY ,StartLine+AccY);      
% % %           
% %             img_sense(:,:,1)= sum(conj(SensMap(:,:,:,slc)).*(img_tmp(:,:,1:nCoil,1)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
% %             img_sense(:,:,2)= sum(conj(SensMap(:,:,:,slc)).*(img_tmp(:,:,1:nCoil,2)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
%             img_tmp=Irecon_EPI_Gfree_fs_comb_vcc(:,:,slc);
%             img_sense=repmat(img_tmp,[1 1 2]);
%             
%             kspace_slice_rc=kspace_slice(:,:,1:nCoil,:);
%             
%             Sens= SensMap(:,:,:,slc);
%             
%             num_iter = 1000;
%             tol = 0.1;
%             lambda = 1.0; 
%             winSize = [3,3];    % local k-space window size 
%             x_k=musselsRecon( img_sense,kspace_slice_rc,Sens,num_iter,tol,winSize,lambda);
% %             mosaic(flipud(mean(abs(sq(x_k)), 3).'), 1, 1, 1, strcat('mussels-',num2str(2),'shots'), [0,1000]), setGcf(.5) 
%        end
%     end
 %%   
    
    
    apodization_para=0.2;
    kcomb_vcc=mrir_fDFT_freqencode(mrir_fDFT_phasencode(Irecon_EPI_Gfree_fs_comb_vcc));
    kcomb=mrir_fDFT_freqencode(mrir_fDFT_phasencode(Irecon_EPI_Gfree_fs_comb));
    
    kapodize_vcc = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb_vcc, 1, apodization_para),  2, apodization_para) ;
    kapodize =mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;
    Irecon_EPI_Gfree_fs_comb_vcc=mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize_vcc));
    Irecon_EPI_Gfree_fs_comb=mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
    
%      mosaic(rsos(Irecon_EPI_Gfree_fs_comb(NRead,:,:),4),1,3,20,['Conv Recon of Mag(dwi) data b',num2str(iReps)]);%caxis([0 50]);
%     mosaic(rsos(Irecon_EPI_Gfree_fs_comb_vcc(NRead,:,:),4),1,3,21,['VCC Recon of Mag(dwi) data b',num2str(iReps)]);%caxis([0 50]);
%      mosaic(rsos(Irecon_EPI_Gfree_fs_comb(:,:,:),4),1,3,20,['Conv Recon of Mag(dwi) data b',num2str(iReps)]);caxis([0 500]);
%     mosaic(rsos(Irecon_EPI_Gfree_fs_comb_vcc(:,:,:),4),1,3,21,['VCC Recon of Mag(dwi) data b',num2str(iReps)]);caxis([0 500]);
%    
    tmp=[];
    if strcmp(TYPE_SAVE, 'SAVE')
       mysave(strcat(savePath,'CPLX_VCC&Conv_Recon_dir',num2str(iReps),'.mat'),Irecon_EPI_Gfree_fs_comb_vcc,Irecon_EPI_Gfree_fs_comb);  
    end 

end
matlabpool close
    clearvars -except iReps
end
  