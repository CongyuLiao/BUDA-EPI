function data_smoothed_i = mrir_image_smooth_2d(data_i, fwhm_voxel)
%MRIR_IMAGE_SMOOTH_2D
%
% img_smoothed = mrir_image_smooth_2d(img, fwhm_voxel)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/07
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  sigma_voxel = fwhm_voxel / ( sqrt(2*log(2))*2 );

  range_voxel = ceil(2 * sigma_voxel);


%  gamma = fwhm_voxel / 2;
%
%  [x_support, y_support] = meshgrid(-range_voxel:+range_voxel);
%  r_support = sqrt(x_support.^2 + y_support.^2);
%  cauchy = (1/pi) * gamma ./ (r_support.^2 + gamma^2);
%
%  cauchy = cauchy / sum(cauchy(:));
%
%  data_smoothed_i = imfilter(data_i, cauchy, 'replicate');


  dims = size(data_i);
  dims(end+1:16) = 1;
  
  
  data_i_2d = reshape(data_i, dims(1), dims(2), []);

  for ind = 1:prod(dims(3:end)),
  
    data_re_smoothed_i_2d = conv2(real(data_i_2d(:,:,ind)), fspecial('gaussian', range_voxel*2, sigma_voxel), 'same');
    data_im_smoothed_i_2d = conv2(imag(data_i_2d(:,:,ind)), fspecial('gaussian', range_voxel*2, sigma_voxel), 'same');
    
    data_smoothed_i_2d(:,:,ind) = complex(data_re_smoothed_i_2d, data_im_smoothed_i_2d);
      
  end;
  
  data_smoothed_i = reshape(data_smoothed_i_2d, dims);
  
  

  return;



  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:




  Ncol = size(data_i, 1);
  Nlin = size(data_i, 2);

  [y, x] = meshgrid([1:Ncol]-Ncol/2,[1:Nlin]-Nlin/2);


  gauss_2d = exp( -(x.*x + y.*y)/(2*sigma_voxel^2) );
  gauss_2d(gauss_2d < eps*max(gauss_2d(:))) = 0;

  area = sum(gauss_2d(:));
  if ( area ~= 0 ),
    gauss_2d  = gauss_2d/area;
  end;

  filt_i = gauss_2d;


  %==--------------------------------------------------------------------==%


  data_k = fft2(data_i, Ncol, Nlin);
  filt_k = fft2(filt_i, Ncol, Nlin);

  dims_data = size(data_k);
  dims_filt = size(filt_k);

  dims_redu = ones(1,length(dims_data));
  dims_redu(1:length(dims_filt)) = dims_filt;


  data_smoothed_i = ifft2(data_k .* repmat(filt_k, dims_data./dims_redu));



  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


