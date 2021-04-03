function [K_epi]= caipi_deblur_CLv1(Kdata, prot, evp, AccZ,PhaseShiftBtwSimulSlices)
K_epi = zeros(size(Kdata));
% detect ky lines
tmp = squeeze(Kdata);
msk = abs(tmp(:,:,1,1)).^(1/3) > 0.001; % Mask to set unacquired points to zero.
ky_idx = find(msk(1,:)~=0);


num_slc = size(Kdata,10)*AccZ;
DThickness = zeros(1,num_slc);
for zz = 1:num_slc
    DThickness(zz) = prot.sSliceArray(zz).dThickness;
end

Kdata_raw = Kdata(:,ky_idx,:,:,:,:,:,:,:,:);

for ii = 1:size(Kdata,10)
    start_slc = ii;
    SliceSep = sum(DThickness(start_slc:(start_slc+(num_slc/AccZ)-1)));  
    K_epi_deblur = CaipirinhaDeblur_CLv2(Kdata_raw, prot, evp, PhaseShiftBtwSimulSlices, SliceSep);
    
end

K_epi(:,ky_idx,:,:,:,:,:,:,:,:) = single(K_epi_deblur);

end