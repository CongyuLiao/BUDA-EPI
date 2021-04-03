function [epi_raw_output, linear_fit_coeff_output, G, epi_img] = mrir_epi_GS_TrainData(meas,KernalX,KernalY)

if isfield(meas, 'patrefscan')
    [epi_raw, acs_raw] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);
else
    epi_raw = meas.data;
    acs_raw = meas.patrefscan;
end
meas.data = [];
meas.patrefscan = [];

s = size(epi_raw); 
s2 = size(meas.data_phascor1d);
epi_raw = mrir_iDFT_freqencode(epi_raw);

% Siemens apply linear phase to shift in-plane FOV: this doesnt work well
% for ramp sampling - so remove before regridding and reapply after - this
% need to be done slice by slice as the shift might be different for diff slice

NormalVec = [meas.prot.sSliceArray(1).sNormal_dSag, meas.prot.sSliceArray(1).sNormal_dCor, meas.prot.sSliceArray(1).sNormal_dTra];
orientation = find(NormalVec == max(NormalVec));
if orientation == 1 % Sagital
    shiftread = [meas.prot.sSliceArray(1:end).sPosition_dSag]; %dCor %dTra
elseif orientation == 2 % Coronal
    shiftread = [meas.prot.sSliceArray(1:end).sPosition_dCor]; %dCor %dTra
else % Transversal
    shiftread = [meas.prot.sSliceArray(1:end).sPosition_dTra]; %dCor %dTra
end

if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    shiftread = reshape(   mrir_image_slice_interleave(reshape(shiftread,1,1,1,1,1,1,1,1,1,[]))  , 1,[]);
end

FOVread = meas.prot.sSliceArray(1).dReadoutFOV *2;% *2 as oversample = 2

%gridding
PhaseUnwindVec = repmat(reshape( single( exp(i*2*pi* (-s(1)/2:(s(1)/2-1)).' *shiftread  /FOVread)), [s(1),1,1,1,1,1,1,1,1,s(10)]), [1, s(2:9), 1]);
epi_raw  = mrir_regrid_trapezoid(epi_raw.*PhaseUnwindVec, meas.prot);
epi_raw = epi_raw.*conj(PhaseUnwindVec);
epi_raw_output= mrir_fDFT_freqencode(sum(epi_raw, 8)); %compress partition
clear PhaseUnwindVec 

PhaseUnwindVecPh = repmat(reshape( single( exp(i*2*pi* (-s2(1)/2:(s2(1)/2-1)).' *shiftread  /FOVread)), [s2(1),1,1,1,1,1,1,1,1,s2(10)]), [1, s2(2:9), 1]);
meas.data_phascor1d = mrir_regrid_trapezoid(meas.data_phascor1d.*PhaseUnwindVecPh, meas.prot);
meas.data_phascor1d = meas.data_phascor1d.*conj(PhaseUnwindVecPh);
clear PhaseUnwindVecPh

%ghost
linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d);
linear_fit_coeff_output = mean(reshape(linear_fit_coeff,2,[],size(meas.data_phascor1d,7)),3);

epi_raw = mrir_artifact_ghost_correct(epi_raw, linear_fit_coeff);
epi_raw = sum(epi_raw, 8); %compress partition

epi_raw = mrir_fDFT_freqencode(epi_raw);

if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_raw = mrir_image_slice_deinterleave(epi_raw);
    epi_raw_output = mrir_image_slice_deinterleave(epi_raw_output);
    linear_fit_coeff_output = reshape(mrir_image_slice_deinterleave(reshape(linear_fit_coeff,[2,1,s(3),1,1,1,1,1,1,s(10)])),2,[]);
end

% original way of doing it. 
% linear_fit_coeff2 = mrir_artifact_ghost_compute(meas.data_phascor1d);
% epi_raw_corr = mrir_artifact_ghost_correct(epi_raw_roft, linear_fit_coeff2);
% epi_raw_grid = mrir_regrid_trapezoid(epi_raw_corr.*PhaseUnwindVec, meas.prot);
% epi_raw_grid = epi_raw_grid.*conj(PhaseUnwindVec);
% epi_2 = mrir_fDFT_freqencode(sum(epi_raw_grid, 8));
% epi_2 = mrir_image_slice_deinterleave(epi_2);

if (~isempty(acs_raw)) % if accelerated acquisition
    epi_acs = mrir_iDFT_freqencode(acs_raw);
    s = size(acs_raw); 
    s2 = size(meas.patrefscan_phascor);
    %gridding
    PhaseUnwindVecACS = repmat(reshape( single( exp(i*2*pi* (-s(1)/2:(s(1)/2-1)).' *shiftread  /FOVread)), [s(1),1,1,1,1,1,1,1,1,s(10)]), [1, s(2:9), 1]);
    epi_acs = mrir_regrid_trapezoid(epi_acs.*PhaseUnwindVecACS , meas.prot);
    epi_acs = epi_acs.*conj(PhaseUnwindVecACS);
    clear PhaseUnwindVecACS

    PhaseUnwindVecACSph = repmat(reshape( single( exp(i*2*pi* (-s2(1)/2:(s2(1)/2-1)).' *shiftread  /FOVread)), [s2(1),1,1,1,1,1,1,1,1,s2(10)]), [1, s2(2:9), 1]);
    meas.patrefscan_phascor = mrir_regrid_trapezoid(meas.patrefscan_phascor.*PhaseUnwindVecACSph, meas.prot);
    meas.patrefscan_phascor = meas.patrefscan_phascor.*conj(PhaseUnwindVecACSph);
    
    %ghost
    linear_fit_coeff_acs = mrir_artifact_ghost_compute(meas.patrefscan_phascor);
    epi_acs = mrir_artifact_ghost_correct(epi_acs, linear_fit_coeff_acs);
    epi_acs = sum(epi_acs, 8); %compress partition
    
    epi_acs = mrir_fDFT_freqencode(epi_acs);
    
    if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
        epi_acs = mrir_image_slice_deinterleave(epi_acs);
    end;
    
    [epi_img, G] = mrir_epi_GRAPPA(epi_raw, epi_acs, meas.prot, meas.evp, KernalX, KernalY, 1, 1);
    
else
    epi_img = mrir_image_crop(mrir_conventional_2d(epi_raw));
    G = [];
end


