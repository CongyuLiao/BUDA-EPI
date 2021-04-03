function img_combine_ada = mrir_array_combine_adaptive(img_uncombined, varargin)
%MRIR_ARRAY_COMBINE_ADAPTIVE
%
% img_combine_ada = mrir_array_combine_adaptive(img_uncombined, noisecovmtx)

% references:
%
%  Walsh DO, Gmitro AF, Marcellin MW.
%  Adaptive reconstruction of phased array MR imagery.
%  Magn Reson Med. 2000 May;43(5):682-690. PMID: 10800033.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/25
% $Id: mrir_array_combine_adaptive.m,v 1.1 2008/01/27 02:51:17 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(img_uncombined, 1)
  Nlin = size(img_uncombined, 2);
  Ncha = size(img_uncombined, 3);

  dims = size(img_uncombined);

  if ( nargin < 2 || isempty(varargin{1}) ),
    Rn = eye(Ncha);
  elseif ( (ndims(varargin{1}) ~= 2) || ...
	   (size(varargin{1}, 1) ~= size(varargin{1}, 2)) ),
    Rn = mrir_array_stats_matrix(Rn, 'cov', 1);
  else,
    Rn = varargin{1};
  end;

  if ( nargin < 3 || isempty(varargin{2}) ),
    Sroi = [1, 1, 1];
  else,
    Sroi = varargin{2};
  end;


  %==--------------------------------------------------------------------==%

  dims(3) = 1;

  % preallocate
  img_combine_ada = zeros(dims);

  rep = 1;
  slc = 1;

  opts_eigs.disp = 0;
  opts_eigs.isreal = 'false';

  % TODO: AC is RSS in special case where sliding window ("ROI") is one
  % pixel, but in general user should be able to define window size in terms
  % of Wcol x Wlin x Wpar.

  % two approaches are possible: a sliding window (such that each pixel's
  % weighting depends on the pixels in its neighborhood, i.e., overlapping
  % ROIs), or a downsampled image (where each pixel's weighting is
  % determined from the macro-pixel to which it belongs).

  METHOD_ROI = 0;

  switch METHOD_ROI,
   case 0,
    % convention adopted by Siemens in functor "AdaptiveCombine.cpp"

   case 1,
    % alternative interpretation

  end;

  iRn = inv(Rn);

  for ii = 1:Ncol,
    for jj = 1:Nlin,

      %                           1,  2, 3, 4, 5, 6,   7, 8, 9,   0
      s = squeeze(img_uncombined(ii, jj, :, 1, 1, 1, rep, 1, 1, slc));
      s = double(s);

      Rs = s*s';

      A = iRn * Rs;
      [w, eigval_max] = eigs(A, 1, 'LM', opts_eigs);

      % [v, d] = eig(A);
      % w = v(:,end);

      img_combine_ada(ii, jj, :, 1, 1, 1, rep, 1, 1, slc) = w' * s;

    end;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_adaptive.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
