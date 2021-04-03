%% load data

load tse1_200x200x21

slice_select = 4;

img = tse1(:,:,slice_select);
img = img / max(img(:));

mosaic(img, 1, 1, 1, 'TSE fully-sampled', [0,0.5])


%% create coil images

load sim_8ch_data

N = size(img);                                          % image size

receive = imresize(b1, N);
receive = receive / max(abs(receive(:)));               % coil sensitivities

num_chan = size(receive,3);                             % no of coils

img_coils = repmat(img, [1,1,num_chan]) .* receive;     % noiseless coil images

PSNR = 200;                                             % add noise with Peak SNR = 200

img_coils = img_coils + (max(abs(img_coils(:))) / PSNR) .* (randn([N, num_chan]) + 1i .* randn([N, num_chan]));     % noisy coil images

img_true = sum(img_coils .* conj(receive), 3) ./ (eps + sum(abs(receive).^2,3));

mosaic(img_true, 1, 1, 2, 'fully-sampled noisy image', [0,.5])
mosaic(receive, 2, 4, 3, 'coil sensitivities', [0,.5])
mosaic(img_coils, 2, 4, 4, 'fully-sampled noisy coil images', [0,.25])


%% ________________________________________________________________________
%% ______________________________ SENSE RECON _____________________________
%% ________________________________________________________________________

m2d = zeros(N);
m2d(1:2:end,1:2:end) = 1;

mosaic(m2d, 1, 1, 4, 'sampling mask', [0,1])

% undersampled coil k-sapce data
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
        
Res = reshape(res, N);
rmse_sense = 100 * norm(Res(:)-img_true(:)) / norm(img_true(:));

mosaic(Res, 1, 1, 5, 'sense recon', [0,.5])
mosaic(5*abs(Res-img_true), 1, 1, 6, ['sense error: ', num2str(rmse_sense), ' % RMSE'], [0,1])
        

