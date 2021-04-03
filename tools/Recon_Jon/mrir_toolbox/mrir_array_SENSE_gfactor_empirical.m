function gfactor = mrir_array_SENSE_gfactor_empirical(i_fullFOV, i_reduFOV, R)


% trim first two time points to remove transients
i_fullFOV(:,:,1:2) = [];

% compute temporal average
i_fullFOV_Tavg = mean(i_fullFOV, 3);
i_reduFOV_Tavg = mean(i_reduFOV, 3);

% compute temporal standard deviation
i_fullFOV_Tstd =  std(i_fullFOV, 0, 3);
i_reduFOV_Tstd =  std(i_reduFOV, 0, 3);


% compute SNR 
i_fullFOV_SNR = i_fullFOV_Tavg ./ i_fullFOV_Tstd;
i_reduFOV_SNR = i_reduFOV_Tavg ./ i_reduFOV_Tstd;


% g-factor is ratio of SNRs, taking into account intrinsic signal loss
% through reduction of samples, sqrt(R)
gfactor = i_fullFOV_SNR ./ i_reduFOV_SNR / sqrt(R);
