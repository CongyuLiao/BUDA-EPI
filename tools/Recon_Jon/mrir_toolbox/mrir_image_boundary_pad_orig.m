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
% [img_pad, mask] = mrir_image_boundary_pad(..., pixel_reach)
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
%  pixel_reach:  number of pixels within mask to use for intensity extrapolation
%  pixel_trim:   number of pixels to remove from boundary of mask

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/may/21
% $Id: mrir_image_boundary_pad.m,v 1.1 2008/05/22 01:23:05 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  thresh = 0.2;

  % this parameter should be roughly the same as the filter kernel size
  pixel_pad = 10;
  pixel_reach = 3;
  pixel_trim = 0;

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
    pixel_reach = varargin{3};
  end;


  if ( nargin >= 5 ),
    pixel_trim = varargin{4};
  end;

  % controls size over which holes in mask are closed up
  pixel_closure = 10;


  %==--------------------------------------------------------------------==%
  % generate pixel mask
  % (PARAMS: thresh, pixel_closure)

  Nframe = 2*max([pixel_pad, pixel_trim, pixel_closure]);

  if ( ~isempty(thresh) ),
    % first mask is a simple intensity threshold
    mask_init = abs(img) > thresh;
  
    if ( DEBUG & nargin >= 2 )
      figure('name', mfilename);
      imagesc(mask_close); axis image; colormap(gray);
      title(sprintf('mask, threshold = %4.1f', thresh));
    end;
  
  end;

  % embed mask in frame to give some room for morphological processing
  mask_frame_init = padarray(mask_init, [Nframe,Nframe], 0, 'both');

  % perform closure to make mask one solid, simply-connected region
  mask_frame_close = imclose(mask_frame_init, strel('disk', pixel_closure));

  % remove frame
  mask_close = mask_frame_close((Nframe+1):(end-Nframe), (Nframe+1):(end-Nframe));
  
  
  %==--------------------------------------------------------------------==%
  % generate band of pixels surrounding mask and overlapping boundary
  % (PARAMS: pixel_trim, pixel_pad)

  % trim off small number of pixels from inner mask containing keeper pixels
  mask_frame_inner = imerode( mask_frame_close, strel('disk', pixel_trim));
  mask_inner = mask_frame_inner((Nframe+1):(end-Nframe), (Nframe+1):(end-Nframe));

  % extend mask out beyond image data -- this is the padding
  mask_frame_dilate = imdilate(mask_frame_close, strel('disk', pixel_pad));

  % find interior extent of band of pixels to be processed
  mask_frame_erode  = imerode( mask_frame_close, strel('disk', pixel_pad));
  mask_erode  = mask_frame_erode((Nframe+1):(end-Nframe), (Nframe+1):(end-Nframe));

  % band is dilated pixels but not pixels left after erosion
  mask_frame_boundary = mask_frame_dilate .* (1-mask_frame_erode);


  mask_boundary = mask_frame_boundary((Nframe+1):(end-Nframe), (Nframe+1):(end-Nframe));


  %==--------------------------------------------------------------------==%
  % extrapolate image intensities into padding
  % (PARAMS: pixel_reach)
  
  
  % median filter to get rid of shot noise
  img_re_median = medfilt2(real(img));
  img_re_boundary = img_re_median .* mask_boundary;

  % dilate the intensities iteratively: large strel applied several times
  % avoids stairstepping artifacts
  img_re_dilate = img_re_boundary;
  for iteration = 1:pixel_reach,
    img_re_dilate = imdilate(img_re_dilate, strel('disk', 5));
  end;

  
  % median filter to get rid of shot noise
  img_im_median = medfilt2(imag(img));
  img_im_boundary = img_im_median .* mask_boundary;

  % dilate the intensities iteratively: large strel applied several times
  % avoids stairstepping artifacts
  img_im_dilate = img_im_boundary;
  for iteration = 1:pixel_reach,
    img_im_dilate = imdilate(img_im_dilate, strel('disk', 5));
  end;

  img_dilate = complex(img_re_dilate, img_im_dilate);
  
  
  % padded image contains original pixels in "mask_inner" and new pixels in
  % "mask_boundary but not mask_inner"
  img_pad = img .* mask_inner  +  img_dilate .* (mask_boundary .* ~mask_inner);


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_boundary_pad.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

