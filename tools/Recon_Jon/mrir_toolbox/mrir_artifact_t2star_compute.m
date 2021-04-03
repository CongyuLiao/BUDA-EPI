function t2star_filter = mrir_artifact_t2star_compute(epi_phascor, MrProt)
%MRIR_ARTIFACT_T2STAR_COMPUTE
%
% t2star_filter = mrir_artifact_t2star_compute(epi_phascor, MrProt)

% note: basic scheme taken from "tdr_fidt2star.m", developed by Doug Greve.
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/16
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  
  if ( size(epi_phascor,1) < 4 ),
    warning('small number (<4) of reference lines may affect filter estimation accuracy');
  end;
  
  % there are more negative k-space indices than positive, so k=0 is at N/2+1.
  k_space_center_index = size(epi_phascor,1)/2 + 1;

  % by default, PE-reversed lines are even, so fit to odd lines only
  fid = epi_phascor(k_space_center_index, 1:2:end);

  Nsamples_to_fit = 8;  % doug uses 5
  
  
  
  X = [ones(nfit,1) -tfid(1:nfit)];
  beta = (inv(X'*X)*X')*log(fid(1:nfit,:));
  T2s = abs(1./beta(2,:));


  

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
