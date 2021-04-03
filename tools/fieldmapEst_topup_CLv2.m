function fieldmapEst_topup_CLv2(save_path,ii_dif,ii_rf)



config_path = './tools/';

fpB0Topup = [save_path, 'img_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'];
% genNii( abs( img_recon_pad ), voxel_size, fpB0Topup)

fpAcqp = [config_path, 'text_ap_pa.txt'];
fpOut = [save_path, 'out_ap_pa_dif',num2str(ii_dif),'rf',num2str(ii_rf)];
fpIout = [save_path, 'epi_unwarp_dif',num2str(ii_dif),'rf',num2str(ii_rf)];
fpField = [save_path, 'fieldmap_dif',num2str(ii_dif),'rf',num2str(ii_rf)];

cmd = ['topup --imain=' fpB0Topup ' --config=', config_path, 'b02b0_topup.cnf --datain=' fpAcqp ' --out=' fpOut ' --iout=' fpIout ' --fout=' fpField]

[status, result] = system(cmd, '-echo');


% load topup results (b0)
% system(['gunzip ', fpIout, '.nii.gz -f'])

% % load field map -> needs to flipped as well
% dt = load_untouch_nii([fpField, '.nii']);
% 
% % img_fieldmap = flipdim(dt.img,1);
% img_fieldmap = dt.img;
end