function img_stat = mrir_image_stats(img, varargin)
%MRIR_IMAGE_STATS
%
% [img_avg, img_std] = mrir_image_stats(img, stattype_str, [S1,S2])
%
%
% stat type can be: 'avg', 'std', 'min', 'max', 'med'


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jun/13
% $Id: mrir_image_stats.m,v 1.1 2008/06/13 17:30:00 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  OPTION_stat = 'avg';

  % use a 3x3 window by default
  [S1,S2] = deal(3);

  if ( nargin >= 2 ),
    OPTION_stat = lower(varargin{1});
  end;

  if ( nargin >= 3 ),

    if ( isempty(varargin{2}) ),

      S1 = size(img, 1);
      S2 = size(img, 2);

    else,

      S1 = varargin{2}(1);
      S2 = varargin{2}(2);

    end;
  end;


  %==--------------------------------------------------------------------==%


  switch OPTION_stat,
   case {'avg', 'std'},

    kernel_avg = fspecial('average', [S1,S2]);
    img_avg = imfilter(img, kernel_avg);

    img_stat = img_avg;

   case {'std'},

    img_demean = img - img_avg;
    img_std = sqrt( imfilter(img_demean.^2, kernel_avg) );

    img_stat = img_std;

   case {'min'},

    img_stat = ordfilt2(img, 1, [S1,S2]);

   case {'max'},

    img_stat = ordfilt2(img, prod(S1,S2), [S1,S2]);

   case {'med'},

    img_stat = medfilt2(img, [S1,S2]);

  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_stats.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
