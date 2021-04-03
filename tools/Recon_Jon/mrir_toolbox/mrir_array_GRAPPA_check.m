function varargout = mrir_array_GRAPPA_check_density(dat, acs)
  

  ind_zeros = find(dat == 0);
  [i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12, i13, i14, i15, i16] = ...
      ind2sub(size(dat), ind_zeros);
  
  if ( ~isempty(i02) ),
    error('sparsity detected in LIN dimension');
  end;
      
  if ( ~isempty(i09) ),
    error('sparsity detected in PAR dimension');
  end;

  
  return;