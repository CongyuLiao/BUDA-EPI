% Step2: estimate field maps using FSL topup
% authors: Congyu Liao, PhD, cyliao@stanford.edu
%          Xiaozhi Cao, PhD, xiaozhic@stanford.edu
% 10/20/2019


clear;clc;close all;


addpath tools
% call_fsl_ubuntu('');

% parallel computing flag:
core_use = 17; %how many cores you wanna use

% save path
save_path = './recon/';

% load protocol and sens map
load([save_path,'prot.mat']);
nRF = 1;  %--1 conventional EPI
          %--5 gSlider factor
nDif =34;

% estimate field maps using FSL topup
 
for ii_rf=1:nRF
    tic
    delete(gcp('nocreate'))
    parpool(core_use) 
    parfor ii_dif=1:nDif   % diffusion directions 1:68
        fieldmapEst_topup_CLv2 (save_path,ii_dif,ii_rf);
        
    end  
    delete(gcp('nocreate'))
    toc
end
