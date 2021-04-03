function mos = mrir_image_mosaic(data, varargin)
%MRIR_IMAGE_MOSAIC  
%
% mos = mrir_image_mosaic(vol)
% mos = mrir_image_mosaic(vol, [row, col])
%
% (wrapper around "vol2mos")

  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/05
% $Id: mrir_image_mosaic.m,v 1.1 2007/02/06 02:21:10 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%  
  
  vol = squeeze(data);
  
  Nslices = size(vol, 3);
  
  
  tilesize = defmossize(Nslices);
  tilelength = max(tilesize);
  
  if ( nargin >= 2 ),
    tile = varargin{1};
  else,
    tile = [tilelength, tilelength];
  end;
  
  
  for ind = 1:size(vol(:,:,:,:), 4),
    mos(:,:,ind) = vol2mos(vol(:,:,:,ind), tile);
  end;
  
  
  D = size(vol);
 
  mos = reshape(mos, [D(1:2).*tile, 1, D(4:end)]);
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/EPI_3D/MATLAB/mrir_image_mosaic.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
