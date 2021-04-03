function [fig, ax] = mrir_display_gfactor(g, varargin);
%MRIR_DISPLAY_GFACTOR  
%
% [fig, ax] = mrir_display_gfactor(g, varargin);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/12
% $Id: mrir_display_gfactor.m,v 1.1 2008/05/31 00:13:58 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  fig = figure('name', inputname(1)); ax = axis; box on;
  imagesc(100./g); axis image; colorbar; caxis([0,100]);
  %  title('SNR loss (%)');


  
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_display_gfactor.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
