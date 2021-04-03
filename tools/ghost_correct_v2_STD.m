function [PAT_ref_PATnavCor_img , PAT_ref_PATnavCor_kspace] = ghost_correct_v2_STD(prot, PAT_ref,PAT_nav, AccY,start_line)

%% ________________________________________________________________________
%% Ghost correction 
%% ________________________________________________________________________

  
disp('Step:  Ghost Correction start')

sz_patrefscan = size(PAT_ref);
nPhases = sz_patrefscan(2);
sy = start_line; % start_line: 3

Ref=zeros(size(PAT_ref));
Ref_RefnavCor=zeros(size(PAT_ref));


Ref(:,sy:AccY:nPhases,:,1,1,1,1,:,1,:)=PAT_ref(:,sy:AccY:nPhases,:,1,1,1,1,:,1,:);
Ref_nav=PAT_nav(:,:,:,:,:,:,:,:,:,:);
lin_fit_RefnavCor = mrir_artifact_ghost_compute_CorrelationMethod(Ref_nav);
Ref_RefnavCor(:,sy:AccY:nPhases,:,1,1,1,1,:,1,:)= mrir_artifact_ghost_correct_CorrelationMethod_v3(Ref(:,sy:AccY:nPhases,:,1,1,1,1,:,1,:), lin_fit_RefnavCor);
clear lin_fit_RefnavCor Ref_nav

data_hybrid = mrir_iDFT_freqencode(Ref_RefnavCor); 

if prot.alRegridMode==1    % if Regrid Mode =1, that means there is no trapezoid gridding.
    data_hybrid_tp=data_hybrid;
else
    data_hybrid_tp = mrir_regrid_trapezoid(data_hybrid, prot);
end
Ref_RefnavCor_kspace = mrir_fDFT_freqencode(data_hybrid_tp);
tmp=sum(Ref_RefnavCor_kspace,8);
PAT_ref_PATnavCor_img = mrir_iDFT_freqencode(mrir_iDFT_phasencode(tmp));
PAT_ref_PATnavCor_kspace=tmp;
disp('Step:  Ghost Correction complete')
end