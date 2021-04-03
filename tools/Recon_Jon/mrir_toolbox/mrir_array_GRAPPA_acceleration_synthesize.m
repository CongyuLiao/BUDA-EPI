function [k_reduFOV, k_ACS] = mrir_array_GRAPPA_acceleration_synthesize(i_fullFOV, R, RawLin, NFirstLin, NRefLin, NFirstRefLin)
%
%
% [k_redu, k_ACS] = mrir_array_GRAPPA_acceleration_synthesize(img, R, RawLin, NFirstLin, NRefLin, NFirstRefLin)
    
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/16
% $Id: mrir_array_GRAPPA_acceleration_synthesize.m,v 1.2 2008/12/16 05:26:23 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( length(R) == 1 ),

    k_fullFOV = mrir_fDFT_freqencode(mrir_fDFT_phasencode(i_fullFOV));
    
    index_datlines = ( NFirstLin : R : (NFirstLin + R * RawLin - 1) ) + 1;
    index_acslines = ( NFirstRefLin : (NFirstRefLin + NRefLin - 1) ) + 1;
    
    k_reduFOV = k_fullFOV(:,index_datlines,:);
    k_ACS     = k_fullFOV(:,index_acslines,:);

  else,
    
    dims = size(i_fullFOV);
    
    % place last dimension (slice or partition) in position 9 following ICE convention
    i_fullFOV_par = reshape(test, [dims(1), dims(2), dims(3), 1, 1, 1, 1, 1, dims(4)]);
    
    k_fullFOV = mrir_fDFT_freqencode(mrir_fDFT_phasencode(mrir_fDFT_phasencode(i_fullFOV_par, 'par'), 'lin'));
    
    index_datlines1 = ( NFirstLin(1) : R(1) : (NFirstLin(1) + R(1) * RawLin(1) - 1) ) + 1;
    index_datlines2 = ( NFirstLin(2) : R(2) : (NFirstLin(2) + R(2) * RawLin(2) - 1) ) + 1;
    index_acslines1 = ( NFirstRefLin(1) : (NFirstRefLin(1) + NRefLin(1) - 1) ) + 1;
    index_acslines2 = ( NFirstRefLin(2) : (NFirstRefLin(2) + NRefLin(2) - 1) ) + 1;
    
    k_reduFOV = k_fullFOV(:,index_datlines,:);
    k_ACS     = k_fullFOV(:,index_acslines,:);
  
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_acceleration_synthesize.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
