function PhaseUnwindVec = mrir_epi_RampReadoutshiftCorr(prot,MatSize)

% Kawin Setsompop 
% Oct 22 2009

% Siemens apply linear phase to shift in-plane FOV: this doesnt work well
% for ramp sampling - so remove before regridding and reapply after - this
% need to be done slice by slice as the shift might be different for diff slice

NormalVec = [prot.sSliceArray(1).sNormal_dSag, prot.sSliceArray(1).sNormal_dCor, prot.sSliceArray(1).sNormal_dTra];
orientation = find(NormalVec == max(NormalVec));
if orientation == 1 % Sagital
    keyboard
    shiftread = [prot.sSliceArray(1:end).sPosition_dSag]; %dCor %dTra dSag
elseif orientation == 2 % Coronal
    keyboard %checksign and 
    shiftread = [prot.sSliceArray(1:end).sPosition_dCor]; %dCor %dTra dSag
else % Transversal: assume read in the Sagital
    shiftread = [prot.sSliceArray(1:end).sPosition_dSag]; %dCor %dTra dSag
end

if ( strcmp(prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    shiftread = reshape(   mrir_image_slice_interleave(reshape(shiftread,1,1,1,1,1,1,1,1,1,[]))  , 1,[]);
end

FOVread = prot.sSliceArray(1).dReadoutFOV *2;% *2 as oversample = 2
PhaseUnwindVec = repmat(reshape( single( exp(i*2*pi* (-MatSize(1)/2:(MatSize(1)/2-1)).' *shiftread/FOVread)), [MatSize(1),1,1,1,1,1,1,1,1,MatSize(10)]), [1, MatSize(2:9), 1]);