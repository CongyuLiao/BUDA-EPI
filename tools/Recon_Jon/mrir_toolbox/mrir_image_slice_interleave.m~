function img_interleaved = mrir_image_slice_interleave(img_deinterleaved)

%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Nslc = mrir_ice_dimensions(img_deinterleaved, 'slc');
  %                             1 2 3 4 5 6 7 8 9        1 2 3 4 5 6
  img_interleaved = img_deinterleaved(:,:,:,:,:,:,:,:,:, ...
                                tdr_sliceorder(Nslc, 1), :,:,:,:,:,:);

  
  return;
  
  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_slice_deinterleave.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
