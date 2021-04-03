function [ res, tflag ] = apply_sense( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels
    
    kspace_coils = reshape(in, [params.N, params.num_chan]);    

    img_coils = zeros([params.N, params.num_chan]);
    for c = 1:params.num_chan
        img_coils(:,:,c) =  ifft2( params.m2d .* kspace_coils(:,:,c) ) * sqrt(prod(params.N));     
    end
    
    Res = sum(img_coils .* conj(params.sens), 3);
    
    res = Res(:);
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
    img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan]) .* params.sens;

    kspace_coils = zeros([params.N, params.num_chan]);
    for c = 1:params.num_chan
        kspace_coils(:,:,c) = params.m2d .* fft2(img_coils(:,:,c)) / sqrt(prod(params.N)); 
    end
    
    res = kspace_coils(:);
    
end

end
