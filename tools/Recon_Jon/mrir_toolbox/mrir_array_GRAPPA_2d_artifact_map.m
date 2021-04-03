function [i_artifact, k_estimate] = mrir_array_GRAPPA_2d_artifact_map(k_fullFOV, G, R1, R2, NLinMeas, NParMeas, varargin)
%MRIR_ARRAY_GRAPPA_2D_ARTIFACT_MAP
%
% mrir_array_GRAPPA_2d_artifact_map(k_fullFOV, G, R1, R2, NLinMeas, NParMeas, varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/10
% $Id: mrir_array_GRAPPA_2d_artifact_map.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  NFirstLin = 1;
  NFirstPar = 1;


  % data collection parameters:

  if ( nargin >= 7 ),
    NFirstLin = varargin{7-6};
  end;

  if ( nargin >= 8 ),
    NFirstPar = varargin{8-6};
  end;


  %==--------------------------------------------------------------------==%

  k_estimate = mrir_array_GRAPPA_2d_artifact(k_fullFOV, G, R1, R2, NLinMeas, NParMeas, NFirstLin, NFirstPar);

  k_artifact = k_fullFOV - k_estimate;
  
  i_artifact = mrir_conventional(k_artifact);



  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_artifact_map.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
