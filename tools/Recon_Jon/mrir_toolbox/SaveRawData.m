function SaveRawData(meas,filename)

if (~isempty(meas.data))
    fname1 =  [filename(1:end-4) '_RawData.data'];
    SaveToDat(meas.data, fname1);
end
if (~isempty(meas.data_phascor1d))
    fname2 =  [filename(1:end-4) '_RawData_nav.data'];
    SaveToDat(meas.data_phascor1d, fname2);
end
meas.data = [];
meas.data_phascor1d = [];

if isfield(meas, 'patrefscan')
    if (~isempty(meas.patrefscan))
        fname3 =  [filename(1:end-4) '_RawACS.data'];
        SaveToDat(meas.patrefscan, fname3);
    end
    if (~isempty(meas.patrefscan_phascor))
        fname4 =  [filename(1:end-4) '_RawACS_nav.data'];
        SaveToDat(meas.patrefscan_phascor, fname4);
    end
    meas.patrefscan = [];
    meas.patrefscan_phascor = [];
end

if isfield(meas, 'smsrefscan')
    if (~isempty(meas.smsrefscan))
        fname5 =  [filename(1:end-4) '_RawSMSref.data'];
        SaveToDat(meas.smsrefscan, fname5);
    end
    if (~isempty(meas.smsrefscan_phascor))
        fname6 =  [filename(1:end-4) '_RawSMSref_nav.data'];
        SaveToDat(meas.smsrefscan_phascor, fname6);
    end
end
meas.smsrefscan = [];
meas.smsrefscan_phascor = [];

if isfield(meas, 'noiseadjscan')   
    if (~isempty(meas.noiseadjscan))
        fname7 =  [filename(1:end-4) '_noiseadjscan.data'];
        SaveToDat(meas.noiseadjscan, fname7);
    end
end

eval(['save ' filename(1:end-4) '_Raw.mat meas'])