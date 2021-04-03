function raw_corr = mrir_artifact_ghost_correct(raw_roft, linear_fit_coeff)
%MRIR_ARTIFACT_GHOST_CORRECT
%
% raw_corr = mrir_artifact_ghost_correct(raw_roft, phase_correct_odd)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/03
% $Id: mrir_artifact_ghost_correct.m,v 1.2 2007/11/14 03:05:16 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

%  if ( mrir_ice_dimensions(raw_roft, 'seg') < 2 ),
%    error('uncorrected data contains only one segment');
%  end;
  
  
  % some sequences (e.g., "ge_functionals") collects reference lines for
  % phase correction every slice.


  % each data dimension will require its own unique phase correction lines,
  % except for the segments dimension (unless truly multi-shot or segmented
  % acquisition) and partitions dimension.

  %   dimension 08: Segments
  %   dimension 09: Partitions

  dims_raw = size(raw_roft);
  dims_raw(end+1:16) = 1;

  if (  ( mrir_ice_dimensions(raw_roft, 'rep') > 1 )  &&  ( prod(dims_raw([3:7,10:16])) ~= size(linear_fit_coeff,2) )  ),
    
    %%% SPECIAL CASE: repetitions
    linear_fit_coeff = reshape(linear_fit_coeff, [2 1 dims_raw(3:6) 1 1 1 dims_raw(10:16)]);
    %                                            1 2 3 4 5 6                                    7 8 9 0 1 2 3 4 5 6
    linear_fit_coeff = repmat(linear_fit_coeff, [1 1 1 1 1 1 mrir_ice_dimensions(raw_roft, 'rep') 1 1 1 1 1 1 1 1 1]);
    linear_fit_coeff = reshape(linear_fit_coeff, 2, []);
    
  end;
   
  
  % if user supplies a single slope and offset for entire data set,
  % replicate the coefficients over all dimensions except Seg and Par
  if ( size(linear_fit_coeff, 2) == 1 );
    linear_fit_coeff = repmat(linear_fit_coeff, [1, prod(dims_raw([3:7,10:16]))]);
  end;

  % values along x-axis over which we evaluate the linear correction
  x_correct = [ones(dims_raw(1),1), ([0:dims_raw(1)-1]-dims_raw(1)/2)'];

  % [ Ncol x prod( N(3:7) N(10:16) ) ]
  phase_correct_odd = x_correct * linear_fit_coeff;

  complex_correct_odd = exp(sqrt(-1) * phase_correct_odd);

  % [ Ncol x 1 x prod( N(3:7) N(10:16) ) ]
  complex_correct_odd = reshape(complex_correct_odd, ...
                                [dims_raw(1), 1, dims_raw(3:7) 1 1 dims_raw(10:16)]);


  dims_reflines = size(complex_correct_odd);
  dims_corrections = dims_reflines(3:end);

  Nd_ref = length(dims_reflines);
  Nd_raw = length(dims_raw);

  % extend list of dimensions for reflines with ones until reaches 1x16
  dims_reflines = cat(2, dims_reflines, repmat(1, [1, (Nd_raw - Nd_ref)]));

  if ( any(dims_reflines > dims_raw) ),
    error('dimension mismatch in phase correction lines');
  end;

  % since various EPI pulse sequence implementations collect phase
  % correction reference lines for each slice/partition, but others might
  % not, we must calculate which of the dimensions must the phase correction
  % reference lines be duplicated. one way to do this is to compare the
  % dimensions of the data volume and the phascor volume. those dimensions
  % for which the data volume contains more data are the ones for which the
  % phase correction reference lines must be duplicated.
  dupl_dims = [1, 1, dims_raw(3:end) - dims_reflines(3:end) + 1];

  % don't duplicate over segments, since we only apply correction to first segment
  dupl_dims(8) = 1;

  % duplicate phase correction lines over those dimensions
  complex_correct_odd_dupl = repmat(complex_correct_odd, dupl_dims);


  % by convention we only calculate the phase discrepancy between the odd
  % and even lines, so we only apply the correction to the odd lines.
  % here we extract the odd lines of the data.
  %                   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  oddlines = raw_roft(:,:,:,:,:,:,:,1,:,:,:,:,:,:,:,:);

  % point-by-point multiplication of odd lines by the phase corrections
  corrected = oddlines .* repmat(complex_correct_odd_dupl, [1,size(oddlines,2)]);

  raw_corr = raw_roft;

  % replace only the odd lines
  %        1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  raw_corr(:,:,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = corrected;


  if ( DEBUG ), mrir_artifact_ghost_correct__PLOT(raw_roft, raw_corr); end;


  return;



%**************************************************************************%
function mrir_artifact_ghost_correct__PLOT(raw_roft, raw_corr)

  
  dims_raw = size(raw_roft);
  
  % visualize first channel
  
  %                                                         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  raw_roft_fDFT = reshape(mrir_fDFT_freqencode(sum(raw_roft(:,:,1,1,1,1,1,:,:,1,1,1,1,1,1,1), 8)), dims_raw(1), []);
  raw_corr_fDFT = reshape(mrir_fDFT_freqencode(sum(raw_corr(:,:,1,1,1,1,1,:,:,1,1,1,1,1,1,1), 8)), dims_raw(1), []);

  figure('Name', mfilename, 'Tag', 'DEBUG_04');
  subplot(2,2,1);
  imagesc(abs(raw_roft_fDFT(:,:)));
  set(gca, 'XTick', 1:size(raw_roft_fDFT,2));
  xlabel('readout line number');
  ylabel('readout sample');
  title('DATA: echo position in k-space');
  subplot(2,2,2);
  imagesc(abs(raw_corr_fDFT(:,:)));
  set(gca, 'XTick', 1:size(raw_corr_fDFT,2));
  xlabel('readout line number');
  ylabel('readout sample');
  title('DATA: echo position corrected');
  subplot(2,2,3);
  imagesc(angle(raw_roft_fDFT(:,:)));
  set(gca, 'XTick', 1:size(raw_roft_fDFT,2));
  xlabel('readout line number');
  ylabel('readout sample');
  title('DATA: echo phase in k-space');
  subplot(2,2,4);
  imagesc(angle(raw_corr_fDFT(:,:)));
  set(gca, 'XTick', 1:size(raw_corr_fDFT,2));
  xlabel('readout line number');
  ylabel('readout sample');
  title('DATA: echo phase corrected');


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_artifact_ghost_correct.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
