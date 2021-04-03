function meas = ReadRawData(filename)

eval(['load ' filename(1:end-4) '_Raw.mat meas'])

fname1 =  [filename(1:end-4) '_RawData.data'];
fname2 =  [filename(1:end-4) '_RawData_nav.data'];
fname3 =  [filename(1:end-4) '_RawACS.data'];
fname4 =  [filename(1:end-4) '_RawACS_nav.data'];
fname5 =  [filename(1:end-4) '_RawSMSref.data'];
fname6 =  [filename(1:end-4) '_RawSMSref_nav.data'];
fname7 =  [filename(1:end-4) '_noiseadjscan.data'];

if ( exist(fname1, 'file') )
    meas.data = ReadFromDat(fname1);
end
if ( exist(fname2, 'file') )
    meas.data_phascor1d = ReadFromDat(fname2);
end
if ( exist(fname3, 'file') )
    meas.patrefscan = ReadFromDat(fname3);
end
if ( exist(fname4, 'file') )
    meas.patrefscan_phascor = ReadFromDat(fname4);
end
if ( exist(fname5, 'file') )
    meas.smsrefscan = ReadFromDat(fname5);
end
if ( exist(fname6, 'file') )
    meas.smsrefscan_phascor = ReadFromDat(fname6);
end
if ( exist(fname7, 'file') )
    dbmeas.noiseadjscan = ReadFromDatNoise(fname7);
end


    
   