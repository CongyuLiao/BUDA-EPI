function F = mrir_array_accelerated_folding_operator_advance(Nlin, R)
%MRIR_ARRAY_ACCELERATED_FOLDING_OPERATOR
%
% F = mrir_array_accelerated_folding_operator_advance(Nlin, R, NlinAliased)

% kawin setsompop June 09
% NlinAliased = number of lines after acceleration
% assumption: NlinAliased = floor(Nlin/R)
%             if Nlin = even, then the centerline (k = 0) is at Nlin/2+1
%             if Nlin = odd, then the centerline (k = 0) is at (Nlin+1)/2
%
% Todo: make it work for non-integer acceleration factor by doing interpolation of k_fullFOV
%**************************************************************************%

  % easy way to calculate folding operator: transform full-FOV image indices
  % to k-space, remove lines skipped during acceleration, then transform
  % back to image space.
  
  NlinAliased = floor(Nlin/R);

  k_fullFOV = mrir_fDFT_freqencode(eye(Nlin));
  
  if rem(NlinAliased,2) == 1 %i.e. odd lines left
      if rem(Nlin,2) == 1 %i.e. odd lines to start with
          Klines = [-floor(NlinAliased/2):floor(NlinAliased/2)]*R  + (Nlin+1)/2;
      else
          Klines = [-floor(NlinAliased/2):floor(NlinAliased/2)]*R  + (Nlin)/2 +1;
      end
  else %i.e. even lines left
     if rem(Nlin,2) == 1 %i.e. odd lines to start with
          Klines = [-(NlinAliased/2):(NlinAliased/2  -1)]*R + (Nlin+1)/2;
      else
          Klines = [-(NlinAliased/2):(NlinAliased/2  -1)]*R + (Nlin)/2 +1;
      end
  end
  
  k_reduFOV = k_fullFOV(Klines, :);
  F = mrir_iDFT_freqencode(k_reduFOV);


