
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% Monday, November 24, 2008 21:43:30 -0500


opt.ReturnStruct = 1;

DATADIR = '/space/pingpong/5/users/jonnyreb/mode_mixing/2008_11_06__mode_mixing_board_test_phantom_HUGO';
meas_noise_normal_BC = read_meas_dat(fullfile(DATADIR, 'meas_MID1757_SNRMap_BC_FID848.dat'), opt);

meas_image_normal_BC = read_meas_dat(fullfile(DATADIR, 'meas_MID1769_gre3D_PD_hugo_BC_FID860.dat'), opt);
meas_image_normal_32 = read_meas_dat(fullfile(DATADIR, 'meas_MID1770_gre3D_PD_hugo_32_FID861.dat'), opt);
img_normal_BC = mrir_iDFT(mrir_iDFT(mrir_iDFT(mean(meas_image_normal_BC.data, mrir_DIM_AVE), 1), 2), 9);
img_BC = mrir_image_crop(img_normal_BC);
mask_thresh = abs(img_BC)>60;
mask = mrir_volume_strip(mask_thresh, 2);
img_32 = mrir_conventional_3d(meas_image_normal_32.data);
surf_BC = isosurface(squeeze(double(mask)), 0.5);



surf_BC = preprocessQ(surf_BC);
surf_BC = preprocessQ(extractpatchCC(surf_BC));
figure; showFace(surf_BC);
redu_BC = reducepatchQ(surf_BC, 0.2);
figure; showFace(redu_BC);
test01 = paint(squeeze(img_32(:,:,01, 1,1,1,1,1, :)), surf_BC, 'nearest');
figure; showVertexValue(surf_BC, test01);
figure; showVertexValue(surf_BC, double(abs(test01)));
camlight
