function [ res, tflag ] = apply_sense_tik( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels

    b = in(1+params.num_chan*prod(params.N):end);    

    kspace_coils = reshape(in(1:params.num_chan*prod(params.N)), [params.N, params.num_chan]);    

    img_coils = zeros([params.N, params.num_chan]);
    for c = 1:params.num_chan
        img_coils(:,:,c) =  ifft2( params.m2d .* kspace_coils(:,:,c) ) * sqrt(prod(params.N));     
    end
    
    Res = sum(img_coils .* conj(params.sens), 3);
    
    res = Res(:) + sqrt(params.lambda) * b;
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
    img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan]) .* params.sens;

    kspace_coils = zeros([params.N, params.num_chan]);
    for c = 1:params.num_chan
        kspace_coils(:,:,c) = params.m2d .* fft2(img_coils(:,:,c)) / sqrt(prod(params.N)); 
    end
    
    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in);
    
end

end
