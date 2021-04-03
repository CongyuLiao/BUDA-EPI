function [meas] = load_BUDA_SMS_data(f_DiffData,iReps,AccZ,evp)

Nreps = 1;
SMSdata = 1;
DataIsVD = 1;

BeginRep = iReps;
[meas] = read_meas_dat_memmap_v3_1_BUDA(f_DiffData,Nreps,BeginRep,SMSdata,0,DataIsVD,AccZ,evp);


% slice index match to gre ref data
if (mod (AccZ,2)==0)
    meas.data = cat(10, meas.data(:,:,:,:,:,:,:,:,:,end),meas.data(:,:,:,:,:,:,:,:,:,1:end-1));
    meas.data_phascor1d = cat(10,meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end),meas.data_phascor1d(:,:,:,:,:,:,:,:,:,1:end-1));
else
    meas.data = meas.data(:,:,:,:,:,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,:);
end



end