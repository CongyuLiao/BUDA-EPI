function PAT_ref_PATnavCor_img  = ghost_correct_pat_ref_v1_FLEET(prot, PAT_ref,PAT_nav, AccY)

%% ________________________________________________________________________
%% Ghost correction for FLEET Ref Data
%% ________________________________________________________________________
%%
%% ________________________________________________________________________
%%

  
disp('Step:  Ghost Correction for FLEET Ref')
TYPE = 'FLEET';
sz_patrefscan = size(PAT_ref);
nPhases = sz_patrefscan(2);
sy = 1; % start line: 3
% AccY=3;
Ref=zeros(size(PAT_ref));
Ref_RefnavCor=zeros(size(PAT_ref));
for R=0: AccY-1
    R
    Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:)=PAT_ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:);
    Ref_nav=PAT_nav(:,R*3+1:R*3+3,:,:,:,:,:,:,:,:);
    lin_fit_RefnavCor = mrir_artifact_ghost_compute_CorrelationMethod(Ref_nav);
    Ref_RefnavCor(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:)= mrir_artifact_ghost_correct_CorrelationMethod_v3(Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:), lin_fit_RefnavCor);
    clear lin_fit_RefnavCor Ref_nav
end
data_hybrid = mrir_iDFT_freqencode(Ref_RefnavCor); 
Ref_RefnavCor=[];
data_hybrid_tp = mrir_regrid_trapezoid(data_hybrid, prot);
Ref_RefnavCor_kspace = mrir_fDFT_freqencode(data_hybrid_tp);
tmp=sum(Ref_RefnavCor_kspace,8);
PAT_ref_PATnavCor_img = mrir_iDFT_freqencode(mrir_iDFT_phasencode(tmp));
disp('Step:  Ghost Correction complete')
end