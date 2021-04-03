function [img_ref_shift, K_ref_shift] = ref_capishift(K_gre,AccZ,PhaseShiftBase)

K_gre = reshape(K_gre,[size(K_gre,1),size(K_gre,2),size(K_gre,3),1,1,1,1,1,1,size(K_gre,4)]);

num_slc = size(K_gre,10);
K_ref_shift = zeros(size(K_gre));
for ii = 1:(num_slc/AccZ)
    start_slc = ii;
    
    K_ref_mslc = K_gre(:,:,:,:,:,:,:,:,:,start_slc:(num_slc/AccZ):num_slc);
    CurrentSliceGroup = [1:(AccZ)];
    K_ref_mslc_shift = zeros(size(K_ref_mslc));
    for ss = 1:AccZ
        K_ref_mslc_shift(:,:,:,:,:,:,:,:,:,ss) = CaipirinhaShift_K_CLv2(K_ref_mslc(:,:,:,:,:,:,1,:,:,ss),CurrentSliceGroup(ss),PhaseShiftBase);
    end
    K_ref_shift(:,:,:,:,:,:,:,:,:,start_slc:(num_slc/AccZ):num_slc) = K_ref_mslc_shift;
end

img_ref_shift = sq(ifft2call(K_ref_shift));
            
            
end