function [img_pad, mask_close] = mrir_image_boundary_pad(img, varargin)
%MRIR_IMAGE_BOUNDARY_PAD
%
% [img_pad, mask] = mrir_image_boundary_pad(img)
%
%
% [img_pad, mask] = mrir_image_boundary_pad(img, mask)
%
% [img_pad, mask] = mrir_image_boundary_pad(img, thresh)
%
% [img_pad, mask] = mrir_image_boundary_pad(..., pixel_pad)
%
% [img_pad, mask] = mrir_image_boundary_pad(..., pixel_trim)
%
%
% options
% -------
%
%  mask:    array size of img to be used as mask
%  thresh:  threshold applied to img to generate simple mask
%
%  pixel_pad:    number of pixels beyond mask to use as pad
%  pixel_trim:   number of pixels to remove from boundary of mask

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/21
% $Id: mrir_image_boundary_pad.m,v 1.3 2008/05/30 23:36:02 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  thresh = 0.2;

  % this parameter should be roughly the same as the filter kernel size
  pixel_pad = 10;

  %
  pixel_trim = 2;


  if ( nargin >= 2 ),
    if ( ~isscalar(varargin{1}) ),
      thresh = [];
      mask_init = varargin{1};
    else,
      if ( nargin >= 2 ),
        thresh = varargin{1};
      end;
    end;
  end;

  if ( nargin >= 3 ),
    pixel_pad = varargin{2};
  end;

  if ( nargin >= 4 ),
    pixel_trim = varargin{3};
  end;


  % controls size over which holes in mask are closed up
  pixel_closure = 10;


  %==--------------------------------------------------------------------==%
  % generate pixel mask
  % (PARAMS: thresh, pixel_closure)

  rss = mrir_array_combine_rss(img);
  
  Nframe = 2*max([pixel_pad, pixel_closure]);

  if ( ~isempty(thresh) ),
    % first mask is a simple intensity threshold
    mask_init = rss > thresh*max(rss(:));

    if ( DEBUG & nargin >= 2 )
      figure('name', mfilename);
      imagesc(mask_init); axis image; colormap(gray);
      title(sprintf('mask init, threshold = %4.1f', thresh));
    end;

  end;

  
  [xx, yy, zz] = meshgrid(-pixel_closure:+pixel_closure);
  sphere_close = sqrt(xx.^2 + yy.^2 + zz.^2) <= pixel_closure;
  
  [xx, yy, zz] = meshgrid(-pixel_trim:+pixel_trim);
  sphere_trim  = sqrt(xx.^2 + yy.^2 + zz.^2) <= pixel_trim;

  
  % embed mask in frame to give some room for morphological processing
  mask_frame_init = padarray(mask_init, Nframe*[1,1,0, 0,0,0,0,0, 1], 0, 'both');
  disp('1');
  
  % perform closure to make mask one solid, simply-connected region
  mask_frame_close = imclose(mask_frame_init, strel('arbitrary', sphere_close));
  disp('2');

  % remove frame
  mask_close = mask_frame_close((Nframe+1):(end-Nframe), (Nframe+1):(end-Nframe));
  disp('3');

  % trim off a layer of pixels from boundary to remove noisy boundary effects
  mask_final = imerode(mask_close, strel('arbitrary', sphere_trim));
  disp('4');
  
  if ( DEBUG )
    figure('name', mfilename);
    imagesc(mask_final); axis image; colormap(gray);
    title(sprintf('mask final, threshold = %4.1f', thresh));
  end;

    
  %==--------------------------------------------------------------------==%
  % generate band of pixels surrounding mask and overlapping boundary
  % (PARAMS: pixel_pad)

  for cha = 1:mrir_ice_dimensions(img, 'cha'),
  
  img_chan = img(:,:,cha);
  img_iter = img_chan;

  mask_grow = mask_final;

  for iteration = 1:pixel_pad,

    % update current image size
    current_size = size(img_iter,1);

    img_re_expand = imresize3(real(img_iter), current_size+2, 'bilinear');
    img_im_expand = imresize3(imag(img_iter), current_size+2, 'bilinear');
    img_expand = complex(img_re_expand, img_im_expand);

    
    img_grow = padarray(img_chan, [1,1,0, 0,0,0,0,0, 1]*iteration, 0, 'both');

    mask_grow = padarray(mask_grow, [1,1,0, 0,0,0,0,0, 1], 0, 'both');

    img_iter = img_grow .* mask_grow  +  img_expand .* (~mask_grow);

    if ( DEBUG ),
      figure('name', mfilename);
      imagesc(abs(img_iter)); axis image; title(sprintf('pad iteration: %3d', iteration));
    end;

  end;

  img_re_pad = imcrop3(real(img_iter), [pixel_pad+1, pixel_pad+1, size(img,1)-1, size(img,2)-1]);
  img_im_pad = imcrop3(imag(img_iter), [pixel_pad+1, pixel_pad+1, size(img,1)-1, size(img,2)-1]);

  img_pad(:,:,cha) = complex(img_re_pad, img_im_pad);
  
  
  if ( DEBUG ),
    figure('name', mfilename);
    imagesc(abs(img_pad(:,:,cha) - img(:,:,cha))); axis image;
    title('expanded region');
  end;
  
  end;
  
  
  return;


%**************************************************************************%
function out = imresize3(img, N, method)
  
  out = imresize(img, [N, 1], method);
  out = ipermute(imresize(permute(out, [2, 1, 3, 4,5,6,7,8, 9]), [N, 1], method), [2, 1, 3, 4,5,6,7,8, 9]);
  out = ipermute(imresize(permute(out, [9, 2, 3, 4,5,6,7,8, 1]), [N, 1], method), [9, 2, 3, 4,5,6,7,8, 1]);


%**************************************************************************%
function out = imcrop3(img, rect)
  
  out = imcrop(img, rect);
  out = ipermute(imcrop(permute(out, [9, 2, 3, 4,5,6,7,8, 1]), [rect(1), 1, rect(3), size(out,2)]), [9, 2, 3, 4,5,6,7,8, 1]);
  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_boundary_pad.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

