function varargout = mrir_array_GRAPPA_check_density(raw)
%MRIR_ARRAY_GRAPPA_CHECK_DENSITY
%
% mrir_array_GRAPPA_check_density(raw)
   
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/nov/13
% $Id: mrir_array_GRAPPA_check_density.m,v 1.2 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  ind_zeros = find(raw == 0);
  [i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12, i13, i14, i15, i16] = ...
      ind2sub(size(raw), ind_zeros);

  if ( ~isempty(i02) && length(i02) > size(raw,1) ),
    warning('sparsity detected in LIN dimension');
  end;

  if ( ~isempty(i09) && length(i09) > size(raw,1) ),
    warning('sparsity detected in PAR dimension');
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_check_density.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

  
  