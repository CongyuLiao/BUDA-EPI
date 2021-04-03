function [epi_raw,linear_fit_coeff] = mrir_epi_GS_RawData(meas)

if isfield(meas, 'patrefscan')
    [epi_raw] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);
    meas.data = []; meas.patrefscan = [];
else
    epi_raw = meas.data;
    meas.data = [];
end

s = size(epi_raw); 

epi_raw = mrir_iDFT_freqencode(epi_raw);

%gridding
epi_raw = mrir_regrid_trapezoid_InplaneShiftCorrect_Interleave(epi_raw, meas.prot);
epi_raw = mrir_fDFT_freqencode(sum(epi_raw, 8)); %compress partition

if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_raw = mrir_image_slice_deinterleave(epi_raw);
end

if nargout == 2
    data_phascor1d = mrir_iDFT_freqencode(meas.data_phascor1d);
    meas.data_phascor1d = [];
    data_phascor1d = mrir_regrid_trapezoid(data_phascor1d, meas.prot);
    
    %ghost
    linear_fit_coeff = mrir_artifact_ghost_compute(mrir_fDFT_freqencode(data_phascor1d));
    clear data_phascor1d 
    
    if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
        linear_fit_coeff = reshape(mrir_image_slice_deinterleave(reshape(linear_fit_coeff,[2,1,s(3),1,1,1,1,1,1,s(10)])),2,[]);
    end
    
end


