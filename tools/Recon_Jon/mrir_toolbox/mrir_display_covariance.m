function [fig, ax] = mrir_display_covariance(covmtx)

  % NOT YET WORKING -- colorbar issue
  

  varvec = diag(covmtx);
  varmtx = diag(varvec);
  stdmtx = sqrt(varmtx);
  rotmtx = inv(stdmtx);
  cofmtx = rotmtx * covmtx * rotmtx;

  cofmtx_dB = 20 * log10(abs(cofmtx));

  
  figure('name', mfilename);
  subplot(1,3,1);
  imagesc(mtx2rgb(abs(covmtx), jet(64)));
  axis image; colorbar;
  title('covariance');

  subplot(1,3,2);
  imagesc(mtx2rgb(abs(cofmtx), arne));
  axis image; colorbar;
  title('correlation coefficient');
  
  subplot(1,3,3);
  imagesc(cofmtx_dB); colormap(hot);
  axis image; colorbar;
  title('correlation coefficient (dB):  20 * log_{10}(cofmtx)');



function mtx_rgb = mtx2rgb(mtx, cmap)
  
  mtx_max = max(mtx(:));
  mtx_min = min(mtx(:));
  
  mtx_scaled = ( mtx - mtx_min ) / ( mtx_max - mtx_min );
  
  mtx_ind = gray2ind(mtx_scaled, size(cmap, 1));
  mtx_rgb = ind2rgb(mtx_ind, cmap);
  
  return;
  

function cmap = arne(varargin)


  cmap = [0            0      0.88547
          0            0            1
          0      0.71922            1
          0      0.88547            1
          0            1            1
          0.71922            1      0.88547
          0.88547            1      0.71922
          1            1            0
          1      0.88547            0
          1      0.71922            0];


  return;
