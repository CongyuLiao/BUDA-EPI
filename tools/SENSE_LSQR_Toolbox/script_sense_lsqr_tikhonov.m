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

num_acs = [24,24];

R = 3;

m2d = rand(N) < (1/R); 

m2d(1+end/2-num_acs(1)/2:end/2+num_acs(1)/2, 1+end/2-num_acs(1)/2:end/2+num_acs(1)/2) = 1;

mosaic(m2d, 1, 1, 4, 'sampling mask', [0,1])

m2d = repmat(m2d, [1,1,num_chan]);


% undersampled coil k-sapce data
kspace_coils = fft2c(img_coils) .* m2d;

% LSQR recon
lsqr_iter = 200;
lsqr_tol = 1e-6;

param.m2d = m2d;
param.sens = receive;
param.N = N;
param.num_chan = num_chan;
param.lambda = 1e-3;        % L2 reg

res = lsqr(@apply_sense_tikc, cat(1, kspace_coils(:), zeros(prod(N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
        
Res = reshape(res, N);
rmse_sense = 100 * norm(Res(:)-img_true(:)) / norm(img_true(:));

mosaic(Res, 1, 1, 5, 'sense recon', [0,.5])
mosaic(5*abs(Res-img_true), 1, 1, 6, ['sense error: ', num2str(rmse_sense), ' % RMSE'], [0,1])
        

