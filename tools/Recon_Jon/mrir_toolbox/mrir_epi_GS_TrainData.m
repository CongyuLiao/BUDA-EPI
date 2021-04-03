function [epi_raw_output, linear_fit_coeff_output, G, epi_img] = mrir_epi_GS_TrainData(meas,KernalX,KernalY)

%grid before  N/2 correction, also apply ramp sampling correction for shift
%inplane

if isfield(meas, 'patrefscan')
    [epi_raw, acs_raw] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);
    meas.data = [];
    meas.patrefscan = [];    
else
    epi_raw = meas.data;
    meas.data = [];
end

s = size(epi_raw); 

epi_raw = mrir_iDFT_freqencode(epi_raw);
data_phascor1d = mrir_iDFT_freqencode(meas.data_phascor1d);
meas.data_phascor1d = [];

%gridding
%epi_raw = mrir_regrid_trapezoid(epi_raw, meas.prot);
epi_raw = mrir_regrid_trapezoid_InplaneShiftCorrect_Interleave(epi_raw, meas.prot);

data_phascor1d = mrir_regrid_trapezoid(data_phascor1d, meas.prot);
%data_phascor1d = mrir_regrid_trapezoid_InplaneShiftCorrect(data_phascor1d, meas.prot); % not needed as effect should be the same for even and odd line

epi_raw_output= mrir_fDFT_freqencode(sum(epi_raw, 8)); %compress partition

%ghost
linear_fit_coeff = mrir_artifact_ghost_compute(mrir_fDFT_freqencode(data_phascor1d));
linear_fit_coeff_output = mean(reshape(linear_fit_coeff,2,[],size(data_phascor1d,7)),3);
clear data_phascor1d

epi_raw = mrir_artifact_ghost_correct(epi_raw, linear_fit_coeff);
epi_raw = sum(epi_raw, 8); %compress partition
epi_raw = mrir_fDFT_freqencode(epi_raw);

if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_raw = mrir_image_slice_deinterleave(epi_raw);
    epi_raw_output = mrir_image_slice_deinterleave(epi_raw_output);
    linear_fit_coeff_output = reshape(mrir_image_slice_deinterleave(reshape(linear_fit_coeff_output,[2,1,s(3),1,1,1,1,1,1,s(10)])),2,[]);
end

if exist('acs_raw','var') % if accelerated acquisition
    s = size(acs_raw);
    s2 = size(meas.patrefscan_phascor);
    epi_acs = mrir_iDFT_freqencode(acs_raw);
    patrefscan_phascor = mrir_iDFT_freqencode(meas.patrefscan_phascor);
    meas.patrefscan_phascor = [];
    
    %gridding
    epi_acs = mrir_regrid_trapezoid_InplaneShiftCorrect(epi_acs, meas.prot);
    %patrefscan_phascor = mrir_regrid_trapezoid_InplaneShiftCorrect(patrefscan_phascor,meas.prot); 
    % not needed as effect should be the same for even and odd line
    patrefscan_phascor = mrir_regrid_trapezoid(patrefscan_phascor, meas.prot);
    
    %ghost
    linear_fit_coeff_acs = mrir_artifact_ghost_compute(mrir_fDFT_freqencode(patrefscan_phascor ));
    clear patrefscan_phascor
    
    epi_acs = mrir_artifact_ghost_correct(epi_acs, linear_fit_coeff_acs);
    epi_acs = sum(epi_acs, 8); %compress partition 
    epi_acs = mrir_fDFT_freqencode(epi_acs);
    
    if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
        epi_acs = mrir_image_slice_deinterleave(epi_acs);
    end;
    
    %[epi_img, G] = mrir_epi_GRAPPA(epi_raw, epi_acs, meas.prot, meas.evp, KernalX, KernalY, 1, 1);
    [epi_img, G] = mrir_epi_GRAPPA(double(epi_raw), double(epi_acs), meas.prot, meas.evp, KernalX, KernalY, 1, 1);
else
    epi_img = mrir_conventional_2d(epi_raw);  
    G = [];

    % debug stuff
    if (0)
        epi_raw = mrir_partial_fourier(mrir_iDFT_freqencode(epi_raw), meas.prot);
        epi = mrir_array_combine(mrir_image_crop(epi_raw),0);
        epi = int16(epi); %o.w. wont save properly for some cases
        save_avw(abs(epi),'img_orderswapPhase','f',[2 2 2 2]);
        unix(['fslview img_orderswapPhase.nii.gz &'])
        save orderswapPhase.mat linear_fit_coeff_output
    end
end


