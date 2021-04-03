function mov = mrir_display_movie(img, dim_row, dim_col, dim_frame)
  %

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/aug/30
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( ndims(img) > 3 ),
  
  dims = size(img);
  
  perm = 1:ndims(img)
  
  ind = find( (perm == dim_row) | (perm == dim_col) | (perm == dim_frame) );
  perm(ind) = [];
  
  perm = [dim_row, dim_col, dim_frame, perm];
  
  img_perm = permute(img, perm);
  img_subset = img( :, :, :, ones(1,ndims(img)-3) );

  else,
    
    img_subset = img;
    
  end;
  
  
  if ( ~isreal(img_subset) ),
    img_subset = abs(img_subset);
  end;
  
  
  for ff = 1:size(img_subset, 3),
    
    indexed = gray2ind(mat2gray((img_subset(:,:,ff))), 64);

    if ( size(indexed,1) < 256 | size(indexed,2) < 256 ),
      
      indexed = imresize(indexed, gray(64), [512, 512], 'Colormap','original');
      
    end;
  
    
    
    mov(ff) = im2frame(indexed, gray(64));
    
  end;
  
  movieview(mov);
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
