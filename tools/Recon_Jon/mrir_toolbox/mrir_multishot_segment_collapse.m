function collapsed = mrir_multishot_segment_collapse(data, varargin)
%MRIR_MULTISHOT_SEGMENT_COLLAPSE
%
% collapsed = mrir_multishot_segment_collapse(data)

% TODO: add functionality to collapse pseudo-segments, e.g., those used
% for phase correction and reflected lines

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/02
% $Id: mrir_multishot_segment_collapse.m,v 1.1 2007/04/06 01:11:39 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  if ( nargin > 1 ),
    evp = varargin{1};
  else,
    evp = [];
  end;

  
  %==--------------------------------------------------------------------==%

  if ( size(data, 8) == 1 ),
    disp('only one segment found; nothing to collapse');
    collapsed = data;
    return;
  end;

  if ( isempty(evp) ),
    % if no evp, then not enough info to perform sparsity checks
    FLAG__SPARSITY_CHECK = 0;
  elseif ( (~isempty(evp.NAFLin) && evp.NAFLin > 1) || ...
	   (~isempty(evp.NAFPar) && evp.NAFPar > 1) ),
    % if data is accelerated, then must skip all (simple) sparsity checks
    FLAG__SPARSITY_CHECK = 0;
  else,
    FLAG__SPARSITY_CHECK = 1;
  end;
  
  
  
  density = nnz(data) / prod(size(data));
  redundancy = 1 / density;

  dimensions = mrir_ice_dimensions(data);

  if ( FLAG__SPARSITY_CHECK && (redundancy ~= dimensions.NSegMeas) ),
    error('sparsity in data does not match number of segments');
  end;

  dims = size(data);
  dims(08) = 1;

  if ( FLAG__SPARSITY_CHECK && (nnz(data) ~= prod(dims)) ),
    error('data can not collapsed along "segments" dimension');
  end;


  % allocate array of NaNs
  collapsed = single(complex(repmat(NaN, dims), repmat(NaN, dims)));

  % initialize collector
  indices_total = [];

  for iseg = 1:dimensions.NSegMeas,

    % extract one segment
    segment = data(:,:,:,:,:,:,:,iseg,:,:,:,:,:,:,:,:);

    % here we assume that none of the data is exactly zero-valued
    indices_segment = find(segment);

    % check to see if any indices in this segment appear in another segment
    if ( any(intersect(indices_total, indices_segment)) ),
      error('repeated data indices detected across segments');
    end;

    % keep a record of all the data indices across segmments
    indices_total = cat(1, indices_total, indices_segment);

    collapsed(indices_segment) = segment(indices_segment);

  end;

  if ( length(unique(indices_total)) < length(indices_total) ),
    error('repetitions detected in collapsed data');
  end;

  if ( FLAG__SPARSITY_CHECK && (any(isnan(collapsed(:)))) ),
    error('some data indices not assigned during collapse');
  end;

  if ( any(isnan(collapsed(:))) ),
    warning('some data indices not assigned during collapse');
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/EPI_3D/MATLAB/mrir_multishot_segment_collapse.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
