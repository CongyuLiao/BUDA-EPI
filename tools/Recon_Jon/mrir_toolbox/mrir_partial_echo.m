function h_full = mrir_partial_echo(k_part, prot, varargin)
%MRIR_PARTIAL_ECHO
%
% h_full = mrir_partial_echo(k_part, prot)

% for use with VIBE data, e.g.:
%
%   meas_MID262_fl3d_vibe_4X4_wOS_192mat_FID18101.dat
%   meas_MID255_fl3d_vibe_128Ch_X12_256mat_FID18094.dat
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/nov/14
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  % default method is zero filling
  METHOD = 'zero-fill';
  if ( nargin >= 3 ),

    % options are 'zero-fill' or 'homodyne'
    METHOD = varargin{1};
  
  end;


  % default omitted k-space data occurs BEFORE echo so TE is shorter in EPI,
  % but sometimes is after echo to shorten total acquisition time in
  % conventional fourier imaging.
  FILL_REGION = 'pre';
  if ( nargin >= 4 ),

    % options are 'pre' or 'post'
    FILL_REGION = varargin{2};
  
    if ( ~ismember(FILL_REGION, {'pre', 'post'}) ),
      error('unrecognized parameter value "%s"; fill region must be either "pre" or "post"', FILL_REGION);
    end;
  
  end;
  
  
  %==--------------------------------------------------------------------==%

  Ncol_full = prot.iNoOfFourierColumns;
  Ncol_part = mrir_ice_dimensions(k_part, 'col');

  if ( Ncol_full == Ncol_part ),

    disp(sprintf('==> [%s]: no partial fourier in COL dimension detected, computing iDFT directly...', mfilename));
    h_full = mrir_iDFT_phasencode(k_part);

    return;
  end;

  k_part_swap_col_lin = permute(k_part, [2,1,3:16]);
  h_full_swap_col_lin = mrir_partial_fourier(k_part_swap_col_lin, struct('lPhaseEncodingLines', Ncol_full), 'zero-fill', 'pre');
  
  h_full = ipermute(h_full_swap_col_lin, [2,1,3:16]);
  

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
