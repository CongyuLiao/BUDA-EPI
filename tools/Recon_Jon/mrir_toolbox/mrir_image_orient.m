function img_orient = mrir_image_orient(img, prot)
%MRIR_IMAGE_ORIENT
%
% img_orient = mrir_image_orient(img, prot)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/19
% $Id: mrir_image_orient.m,v 1.2 2008/04/20 00:02:07 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( sum([prot.sSliceArray(1).sNormal_dSag,
            prot.sSliceArray(1).sNormal_dCor,
            prot.sSliceArray(1).sNormal_dTra]) ~= 1 ),
    error('incorrect slice normal detected');
  end;


  %==--------------------------------------------------------------------==%

  if ( prot.sSliceArray(1).sNormal_dSag ),
    img_orient = mrir_image_orient_SAG(img);
  end;
  if ( prot.sSliceArray(1).sNormal_dCor ),
    img_orient = mrir_image_orient_TRA(img);
  end;
  if ( prot.sSliceArray(1).sNormal_dTra ),
    img_orient = mrir_image_orient_TRA(img);
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_orient.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: