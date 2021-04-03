function [Img_gre, whtcc]= load_GRE_data_prewht(fGRE,N,nc)

dt = read_meas_dat(fGRE);
gre = dt.data;
gre = mrir_image_slice_deinterleave( gre );
% remove OS
gre = ifft2call(permute(squeeze(gre), [1 2 4 3]));
gre = fft2call(gre(1+end/4:3*end/4,:,:,:,:));
[n(1), n(2), ~, ~] = size(gre);
% zero pad patref k-space
rref = ifft2call(padarray(gre, (N-n)/2));


% load noise matrix
dt1 = mapVBVD(fGRE);
nmat = permute(squeeze(dt1{1}.noise()), [2, 1, 3, 4]);
nmat = reshape(nmat, size(nmat, 1), prod(size(nmat))/size(nmat, 1));


[rx, ry, rz, rc] = size(rref);
% SVD to get coil compression matrix.
rref = reshape(rref, rx * ry * rz, rc).';
[u, s, ~] = svd(rref, 'econ');
u = u(:, 1:nc);
s = diag(s);

% Coil compressing noise.
nmat = u' * nmat;

% Estimating whitening matrix.
covm = (nmat * nmat')/(size(nmat, 2) - 1);
whmt = inv(chol(covm, 'lower'));
whmt=whmt./max(max(abs(whmt)));

% Joint coil compression and whitening matrix.
whtcc = whmt * u';

% Whiten and coil compress reference data.
rref = whtcc * rref;
% rref = u'*rref;
rref = rref.';
Img_gre = reshape(rref, rx, ry, rz, nc);
Img_gre = permute( Img_gre, [1 2 4 3]);



end