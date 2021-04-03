function [Img_gre, K_gre]= load_GRE_data(fGRE,N)

dt = read_meas_dat(fGRE);
gre = dt.data;
gre = mrir_image_slice_deinterleave( gre );
% remove OS
gre = ifft2call(sq(gre));
gre = fft2call(gre(1+end/4:3*end/4,:,:,:,:));
[n(1), n(2), ~, ~] = size(gre);
% zero pad patref k-space
Img_gre = ifft2call(padarray(gre, (N-n)/2));
K_gre = padarray(gre, (N-n)/2);


end