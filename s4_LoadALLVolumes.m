% Step4: load whole-brain volumes and do real-valued diffusion processing 
% authors: Congyu Liao, PhD, cyliao@stanford.edu
%          Xiaozhi Cao, PhD, xiaozhic@stanford.edu
% 10/20/2019



% load whole brain data
clear;clc;close all;
addpath(genpath('./tools'))
addpath tools
save_path='./recon/';

nRF = 1;  % --1 conventional EPI
          %--5 gSlider factor
nDif =34;  % 34 diffusion directions

ii = 0;
for ii_dif = 1:nDif
    for ii_rf = 1:nRF
        ii_dif
        ii_rf  
        ii= ii+1;
        load([save_path, 'img_buda_rf', num2str(ii_rf),'_dif', num2str(ii_dif),'.mat'])        
        img_buda_all(:,:,:,ii,:)=permute(img_msl,[1 2 4 3]);
        clear img_sense img_hybrid_sense img_msl
        
        
    end
end
%% do real value diffusion
addpath ('library')
for ii=1:size(img_buda_all,5)
    ii
    img_buda_real(:,:,:,:,ii)=single(RealDiffusion_lowRes(img_buda_all(:,:,:,:,ii),1,1,0));

end
img_buda_real=mean(img_buda_real,5);

disp('saving')
save([save_path, 'buda_real_TR',num2str(nDif*nRF),'.mat'], 'img_buda_real','-V7.3')


% clearvars -except load_path img_hybrid_sense_real img_buda_real img_ap_real img_pa_real  


