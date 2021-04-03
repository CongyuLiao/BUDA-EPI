function ph = mrir_artifact_t2star_plot(epi_phascor, MrProt)
%MRIR_ARTIFACT_T2STAR_PLOT
%
% ph = mrir_artifact_t2star_plot(epi_phascor, MrProt)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/24
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % there are more negative k-space indices than positive, so k=0 is at N/2+1.
  k_space_center_index = size(epi_phascor,1)/2 + 1;

  % by default, PE-reversed lines are even, so fit to odd lines only
  fid = squeeze(epi_phascor(k_space_center_index, 1:2:end,:,:,:,:,:,:,:));

  figure('Name', mfilename);
  ph = plot(abs(fid));
  legend(cellstr(reshape(sprintf('chan %3d', 1:size(fid,2)), [8,size(fid,2)])'));


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
