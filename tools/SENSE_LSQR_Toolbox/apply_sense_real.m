function [ res, tflag ] = apply_sense_tikc( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels


    kspace_re = reshape(in(1:end/2), [params.N, params.num_chan]);    
    kspace_im  = reshape(in(1+end/2:end), [params.N, params.num_chan]);    

    img_coils = ifft2c( kspace_re .* params.m2d );
    Re = real( sum(img_coils .* conj(params.sens), 3) );
    
    img_coils = ifft2c( kspace_im .* params.m2d );
    Im = imag( sum(img_coils .* conj(params.sens), 3) );

    res = Re(:) + Im(:);
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
    img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan]) .* params.sens;

    kspace_coils = fft2c(img_coils) .* params.m2d;
    
    res = cat(1, real(kspace_coils(:)), imag(kspace_coils(:)));
    
end

end
