function [img_RAS, aspect] = mrir_image_vox2ras(img_VOX)
  
  
  nSLI = size(img_RCS, 10);
  nPAR = size(img_RCS, 09);
  
  
  % 2D: SLI dim moves to position 3
  img_RCS = squeeze(permute(img_VOX, [1, 2, 10, 3:9, 11:16]));

  % 3D: PAR dim moves to position 3
  img_RCS = squeeze(permute(img_VOX, [1, 2, 09, 3:8, 10:16]));