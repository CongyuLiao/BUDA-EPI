function mrir_artifact_ghost_compute__DEBUG(phascor, phasediff, linear_fit_coeff, ...
                                            ref_REV, ref_FWD, dims)

  % "DEBUG" must be set if this function is called. to avoid duplicate
  % debugging figures, temporarily disable the debugging features (clearing
  % removes global variable in local scope)
%  clear DEBUG

  x_refline = [ones(size(phasediff,1),1), ...
               ([0:size(phasediff,1)-1]-size(phasediff,1)/2)'];


  x_correct = [ones(dims(1),1), ([0:dims(1)-1]-dims(1)/2)'];
  phase_correct_REV = x_correct * linear_fit_coeff;

  complex_correct_REV = exp(sqrt(-1) * phase_correct_REV);
  complex_correct_REV = reshape(complex_correct_REV, [dims(1), 1, dims(3:end)]);

  projections = mrir_iDFT_freqencode(phascor);
%
%  dims = ones(16,1);
%  dims(8) = 2;
%  projections_seg = repmat(zeros(size(projections)), dims);
%  %               1       2 3 4 5 6 7 8
%  projections_seg(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = projections(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
%  projections_seg(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:) = projections(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
%
  projections_seg = projections;
 
  global DEBUG
  DEBUG = 0;
  corrections = sum(mrir_artifact_ghost_correct(projections_seg, linear_fit_coeff), 8);
  DEBUG = 1;
  
  reflines_correct = mrir_fDFT_freqencode(corrections);

  reflines_recon = mrir_iDFT_phasencode(corrections);
  reflines_recon = mrir_image_crop(reflines_recon);

  reflines_recon = mrir_array_combine_rss(reflines_recon);

  reflines_recon = mrir_intensity_scale_bin(reflines_recon, 4095);
  reflines_recon = permute(flipdim(reflines_recon, 2), [2 1 3:16]);

  Nref = min([8, size(projections,2)]);


  %==--------------------------------------------------------------------==%

%  figure('Name', mfilename);
%  plot(abs((phascor(:,:, 1))));
%  legend(cellstr(reshape(sprintf('line %3d', 1:size(phascor,2)), [8,size(phascor,2)])'));
%  title('echoes in k-space');

%  figure('Name', mfilename); axes; hold on;
%  plot(x_correct(:,2), (angle(ref_REV(:,1))) ,'k');
%  plot(x_correct(:,2), (angle(ref_FWD(:,1))), 'm');
%  plot(x_correct(:,2), (angle(ref_REV(:,1)./ref_FWD(:,1))), 'g');
%  legend('ref_REV', 'ref_FWD', 'diff');
%  title('averaged phase after readout ifft');

  figure('Name', mfilename, 'Tag', 'DEBUG_01');
  subplot(3,3,2); hold on;
  plot(x_correct(:,2), unwrap(angle(ref_REV(:,1))) ,'k');
  plot(x_correct(:,2), unwrap(angle(ref_FWD(:,1))), 'm');
  plot(x_correct(:,2), unwrap(angle(ref_FWD(:,1)./ref_REV(:,1))), 'g');
  lh = legend('ref_REV', 'ref_FWD', 'diff', 'Location', 'NorthEastOutside');
  set(lh, 'Interpreter', 'none')
  title('averaged unwrapped phase after readout ifft');

  subplot(3,3,5); hold on;
  plot(x_refline(:,2), phasediff(:,1), 'ro');
  plot(x_correct(:,2), phase_correct_REV(:,1), '.-'); axis equal; grid on;
  legend('phase diff', 'linear fit', 'Location', 'NorthEastOutside');

  subplot(3,3,8); hold on;
  plot(angle(ref_FWD(:,1)), 'r');
  plot(angle(ref_REV(:,1)), 'b');
  plot(angle(ref_REV(:,1).*complex_correct_REV(:,1)), 'g');
  lh = legend('ref_FWD', 'ref_REV', 'ref_REV corrected', 'Location', 'NorthEastOutside');
  set(lh, 'Interpreter', 'none')
  
  phascor_collapse = sum(phascor, 8);
  
  %trial = min(4, size(phascor,3));
  trial = 1;
  figure('Name', mfilename, 'Tag', 'DEBUG_02');
  subplot(2,2,1);
  imagesc(abs(phascor_collapse(:,:, trial)));
  set(gca, 'XTick', 1:size(phascor,2));
  xlabel('reference line number');
  ylabel('readout sample');
  title('echo position in k space');
  subplot(2,2,2);
  imagesc(abs(reflines_correct(:,:, trial)));
  set(gca, 'XTick', 1:size(reflines_correct,2));
  xlabel('reference line number');
  ylabel('readout sample');
  title('echo position corrected');
  subplot(2,2,3);
  imagesc(angle(phascor_collapse(:,:, trial)));
  set(gca, 'XTick', 1:size(phascor,2));
  xlabel('reference line number');
  ylabel('readout sample');
  title('echo phase in k-space');
  subplot(2,2,4);
  imagesc(angle(reflines_correct(:,:, trial)));
  set(gca, 'XTick', 1:size(reflines_correct,2));
  xlabel('reference line number');
  ylabel('readout sample');
  title('echo phase corrected');

%  figure('Name', mfilename);
%  plot(unwrap(angle(fftshift(reflines_correct(:,1:Nref, trial), 1))));
%  legend(cellstr(reshape(sprintf('line %3d', 1:Nref), [8,Nref])'));
%  axis equal;
%  axis([-5, 15, -2*pi, +2*pi])


  figure('Name', mfilename, 'Tag', 'DEBUG_03');
  imagesc(abs(reflines_recon(:,:,1)));
  colormap(gray);
  set(gca, 'YTick', [1:size(reflines_recon,2)]);
  xlabel('readout direction');
  ylabel('reference lines');
  title('reconstruction of reference lines: projection of object');

  return;
