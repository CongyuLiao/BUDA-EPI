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

tic
    res = lsqr(@apply_sense, kspace_coils(:), lsqr_tol, lsqr_iter, [], [], [], param);  
toc
    
Res = reshape(res, N);
rmse_sense = 100 * norm(Res(:)-img_true(:)) / norm(img_true(:));

mosaic(Res, 1, 1, 5, 'sense recon', [0,.5])
mosaic(5*abs(Res-img_true), 1, 1, 6, ['sense error: ', num2str(rmse_sense), ' % RMSE'], [0,1])
        

%% split data consistency

lambda = eps;
alpha = 1;

M3D = repmat(m2d, [1,1,num_chan]);
I3D = ones([N, num_chan]) * alpha;
I2D = ones(N) * lambda / alpha;
CtC = sum(abs(receive).^2, 3);

zf = zeros([N, num_chan]);
for c = 1:num_chan
    zf(:,:,c) = ifft2(kspace_coils(:,:,c)) * sqrt(prod(N));
end

% x = zeros(N);
x = sum(conj(receive) .* zf, 3) ./ (eps + CtC);

y = zeros([N, num_chan]);

tic
for t = 1:150

    % update y = C * x
    Cx = repmat(x, [1,1,num_chan]) .* receive;
    FCx = zeros([N, num_chan]);
    
    for c = 1:num_chan
        FCx(:,:,c) = fft2(Cx(:,:,c)) / sqrt(prod(N));
    end
    
    temp = (M3D .* kspace_coils + alpha .* FCx) ./ (M3D + I3D);
    
    for c = 1:num_chan
        y(:,:,c) = ifft2(temp(:,:,c)) * sqrt(prod(N));
    end
    
    % update x
    x = sum(conj(receive) .* y, 3) ./ (CtC + I2D);
    
    if mod(t, 20) == 0
        mosaic(x, 1, 1, 10, 'split recon', [0,.5])
    end
end
toc

rmse_split = 100 * norm(x(:)-img_true(:)) / norm(img_true(:));

mosaic(5*abs(x-img_true), 1, 1, 11, ['split error: ', num2str(rmse_split), ' % RMSE'], [0,1])




