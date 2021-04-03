function sensitivity_trunc = mrir_array_modemix_compress_pruessmann(sensitivity, chan)
  
  
  
  sigcormtx = mrir_array_stats_matrix(sensitivity, 'cor');
  
  [U, S, V] = svd(sigcormtx);
  
  
  sensitivity_mixed = mrir_array_coord_transform(sensitivity, U');
  
  %                                     1 2      3 4 5 6 7 8 9 0 1 2 3 4 5 6
  sensitivity_trunc = sensitivity_mixed(:,:,1:chan,:,:,:,:,:,:,:,:,:,:,:,:,:);
  
  
  return;