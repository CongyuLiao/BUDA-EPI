function [ res, tflag ] = apply_sense3D( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels
    
    kspace_coils = reshape(in, [params.N, params.num_chan]);    
    
    
    
    if isfield(params, 'shift')    
    
        if params.shift

            img_coils = fftshift(ifft(fftshift(kspace_coils .* params.m3d, 3), [], 3), 3);
            img_coils = fftshift(ifft2(fftshift(img_coils))) * params.fft_scale;

        else

            img_coils = ifft(kspace_coils .* params.m3d, [], 3);
            img_coils = ifft2(img_coils) * params.fft_scale;     

        end
        
    else
        
        % fftshift by default
        img_coils = fftshift(ifft(fftshift(kspace_coils .* params.m3d, 3), [], 3), 3);
        img_coils = fftshift(ifft2(fftshift(img_coils))) * params.fft_scale;
        
    end
       
        
        
    
    
    Res = sum(img_coils .* params.csens, 4);
    
    res = Res(:);
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
    img_coils = repmat(reshape(in, params.N), [1,1,1,params.num_chan]) .* params.sens;


    if isfield(params, 'shift')    

        if params.shift
            
            kspace_coils = fftshift(fft(fftshift(img_coils, 3), [], 3), 3);
            kspace_coils = fftshift(fft2(fftshift(kspace_coils))) .* params.m3d / params.fft_scale;
            
        else
            
            kspace_coils = fft(img_coils, [], 3);
            kspace_coils = fft2(kspace_coils) .* params.m3d / params.fft_scale;

        end
    
    else
        
        % fftshift by default
        kspace_coils = fftshift(fft(fftshift(img_coils, 3), [], 3), 3);
        kspace_coils = fftshift(fft2(fftshift(kspace_coils))) .* params.m3d / params.fft_scale;

    end
    
    
    res = kspace_coils(:);
    
    fprintf('+')
    
end

end


