function [ res, tflag ] = apply_sense_fixedPhase( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels

    b = in(1 + params.num_chan * prod(params.N) * params.num_shot : end);    

    kspace_coils = reshape(in(1 : params.num_chan * prod(params.N) * params.num_shot), [params.N, params.num_chan, params.num_shot]);    

    img_coils = ifft2c2( kspace_coils .* params.m2d );
        
%     Res = sum(img_coils .* conj(params.sens), 3);
  
    Res = sum(sum(img_coils .* conj(params.phs_sens), 3), 4);
    
    res = Res(:) + sqrt(params.lambda) * b;
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT

    % params.phs_sens : sens .* phs for each shot 
    
    img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan, params.num_shot]) .* params.phs_sens;
            
%     img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan]) .* params.sens;

    kspace_coils = fft2c2(img_coils) .* params.m2d;
    
    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in);
    
end

end
