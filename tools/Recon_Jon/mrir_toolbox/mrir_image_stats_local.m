function img_stat = mrir_image_stats_local(img, varargin)
%MRIR_IMAGE_STATS_LOCAL
%
% [img_avg, img_std] = mrir_image_stats_local(img, [S1, S2])

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jun/13
% $Id: mrir_image_stats_local.m,v 1.1 2008/06/13 17:05:15 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % use a 3x3 window by default
  [S1,S2] = deal(3);

  if ( nargin >= 2 ),

    S1 = varargin{1}(1);
    S2 = varargin{1}(2);

  end;


  %==--------------------------------------------------------------------==%

  kernel_avg = fspecial('average', [S1,S2]);


  img_avg = imfilter(img, kernel_avg);


  img_demean = img - img_avg;

  img_std = sqrt( imfilter(img_demean.^2, kernel_avg) );



  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_stats_local.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
