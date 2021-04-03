function [sensitivity, img_AC, img_BC] = mrir_sensitivity_map_prescan_normalize(meas)
%MRIR_SENSITIVITY_MAP_PRESCAN_NORMALIZE  recon sensitivities from "AdjCoilSens" data
%
% sensitivity = mrir_sensitivity_map_prescan_normalize(meas_struct)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/07
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  raw_AC = meas.data(:,:,:,1,:,:,:,:,:,:);
  raw_BC = meas.data(:,:,1,2,:,:,:,:,:,:);

  img_AC = mrir_conventional_3d(raw_AC);
  img_BC = mrir_conventional_3d(raw_BC);

   
  
  sensitivity = img_AC ./ repmat(img_BC+100, size(raw_AC)./size(raw_BC));
  
  


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
%
%
%"TotalNoOfCoilsInProt"
%34
%"CoilSetup"
%3
%
%asCoilSelectMeas[0].asList[0].sCoilElementID.tCoilID = "32Channel_MGH"
%asCoilSelectMeas[1].asList[0].sCoilElementID.tCoilID = "Body"
%asCoilSelectMeas[1].asList[0].sCoilElementID.tElement = "BC2"
%asCoilSelectMeas[1].asList[1].sCoilElementID.tCoilID = "Body"
%asCoilSelectMeas[1].asList[1].sCoilElementID.tElement = "BC"
