%% load data


 
load /cluster/kawin/Bo/sense_recon
load /cluster/kawin/Bo/sense_map


receive = sense_map;

img_coils = fftshift(ifft2(fftshift(kData)));

% img_coils = padarray(img_coils, [16,0,0]);



N = size(img_coils(:,:,1));                                          % image size
num_chan = size(receive,3);                             % no of coils



img_true = sum(img_coils .* conj(receive), 3) ./ (eps + sum(abs(receive).^2,3));


 

%%  



mosaic(img_true, 1, 1, 2, 'fully-sampled noisy image', [0,5e-3])
mosaic(receive, 2, 4, 3, 'coil sensitivities', [0,.5])
mosaic(img_coils, 2, 4, 4, 'fully-sampled noisy coil images', [0,1e-3])




%% ________________________________________________________________________
%% ______________________________ SENSE RECON _____________________________
%% ________________________________________________________________________


m2d = zeros(N);

m2d(1:2:end,1:1:end) = 1;
% m2d(1+end/2-25:end/2+25,:) = 1;


% m2d(1:1:end,1:2:end) = 1;
% m2d(:,1+end/2-25:end/2+25) = 1;


m2d = fftshift(m2d);


mosaic(m2d, 1, 1, 4, 'sampling mask', [0,1])




% undersampled coil k-space data
kspace_coils = zeros([N, num_chan]);
for c = 1:num_chan
    kspace_coils(:,:,c) = fft2(img_coils(:,:,c)) .* m2d / sqrt(prod(N));
end


% LSQR recon
lsqr_iter = 100;
lsqr_tol = 1e-6;

param.m2d = m2d;
param.sens = receive;
param.N = N;
param.num_chan = num_chan;



res = lsqr(@apply_sense, kspace_coils(:), lsqr_tol, lsqr_iter, [], [], [], param);  
        


Res = reshape(res, N);.
rmse_sense = 100 * norm(Res(:)-img_true(:)) / norm(img_true(:));


mosaic(img_true, 1, 1, 1, 'sense R=1', [0,5e-3])
mosaic(Res, 1, 1, 200, 'sense recon', [0,5e-3])
mosaic(recon_sense / sqrt(prod(N)), 1, 1, 3, 'sense recon', [0,5e-3])


% mosaic(100 * (Res - recon_sense / sqrt(prod(N))), 1, 1, 13, 'sense recon', [0,5e-3])



mosaic(5*abs(Res-img_true), 1, 1, 206, ['sense error: ', num2str(rmse_sense), ' % RMSE'], [0,5e-3])





%%


c = 0.6;
num_acs = 16;

 writecfl('img', permute(single(img_coils), [4,1,2,3]))

 system(['fft 7 ', 'img ', 'kspace'])
 
 system(['ecalib -c ', num2str(c), ' -r ', num2str(num_acs), ' ', 'kspace ', 'calib'])
 
 system(['slice 4 0 ', 'calib ',  'sens'])
    
     
     
     
 receive = squeeze( single(readcfl(['sens'])) );

     
     

        

