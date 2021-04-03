function [img_sense_topup] = applytopup_APPA(img_recon_pad,save_path,voxel_size,ii_dif,ii_rf)
img_sense_pad = img_recon_pad;
config_path = '/autofs/cluster/kawin/Congyu/CODE/gSlider_BUDA/Recon_Data/';

fpAcqp = [config_path, 'acq_param.txt'];
fpOut = [save_path, 'out_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf)];


system(['rm ', save_path, 'my_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])

system(['rm ', save_path, 'my_good_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])


genNii( real(img_sense_pad(:,:,:,1)), voxel_size, [save_path, 'my_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
genNii( imag(img_sense_pad(:,:,:,1)), voxel_size, [save_path, 'my_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])

genNii( real(img_sense_pad(:,:,:,2)), voxel_size, [save_path, 'my_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
genNii( imag(img_sense_pad(:,:,:,2)), voxel_size, [save_path, 'my_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])



cmd = ['applytopup --imain=', save_path, 'my_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ...
    ' --method=jac --out=', save_path, 'my_good_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --interp=spline'];

[status, result] = system(cmd, '-echo');


cmd = ['applytopup --imain=', save_path, 'my_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ...
    ' --method=jac --out=', save_path, 'my_good_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --interp=spline'];

[status, result] = system(cmd, '-echo');



cmd = ['applytopup --imain=', save_path, 'my_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ...
    ' --method=jac --out=', save_path, 'my_good_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --interp=spline'];

[status, result] = system(cmd, '-echo');


cmd = ['applytopup --imain=', save_path, 'my_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ...
    ' --method=jac --out=', save_path, 'my_good_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), ' --interp=spline'];

[status, result] = system(cmd, '-echo');


system(['gunzip ', save_path, '*.gz -f'])


dt_re = load_untouch_nii([save_path, 'my_good_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), '.nii']);
dt_im = load_untouch_nii([save_path, 'my_good_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), '.nii']);

img_1 = dt_re.img + 1i * dt_im.img;



dt_re = load_untouch_nii([save_path, 'my_good_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf), '.nii']);
dt_im = load_untouch_nii([save_path, 'my_good_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf), '.nii']);

img_2 = dt_re.img + 1i * dt_im.img;


img_sense_topup(:,:,:,:) = permute(cat(4, img_1, img_2), [1,2,4,3]);

system(['rm ', save_path, 'my_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])

system(['rm ', save_path, 'my_good_blipup_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipup_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipdown_re_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])
system(['rm ', save_path, 'my_good_blipdown_im_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'])



% mosaic(img_sense_topup(:,:,1,:,1), 6, 8, 1, 'AP topup ', genCaxis(img_sense_topup(:,:,:,:,1)), 90)
% mosaic(img_sense_topup(:,:,2,:,1), 6, 8, 2, 'PA topup ', genCaxis(img_sense_topup(:,:,:,:,1)), 90)
% 
% 



end