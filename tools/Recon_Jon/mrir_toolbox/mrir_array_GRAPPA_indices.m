function [datlines, acslines, reflines, datindex, acsindex] = ...
    mrir_array_GRAPPA_indices(R, datfirst, NDat, acsfirst, NAcs)
%MRIR_ARRAY_GRAPPA_INDICES
%
% [datlines, acslines, reflines, datindex, acsindex] = ...
%    mrir_array_GRAPPA_indices(R, datfirst, NDat, acsfirst, NAcs)
  
% datlines: indices of lines in sparse data array of data
% acslines: indices of lines in sparse data array of ACS measurements
% reflines: indices of lines in sparse data array of extra reference lines
% datindex:
% acsindex:


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/apr/06
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  
  datlines = datfirst : R : NDat;
  datindex(datlines) = 1:length(datlines);
 
  if ( NAcs <= 1 ),
    % if no acceleration in partition direction...
    acslines = 1;
  else,
    if ( acsfirst == 0 ),  % this happens with separate ACS lines
      acsfirst = 1;
    end;
    
    % ASSUMPTION: ACS lines are always contiguous
    %acslines = 1:60;  % cxz_test
    acslines = ( 1:NAcs ) + acsfirst - 1; % cxz_test_original
    acsindex(acslines) = 1:length(acslines);
  end;
  
  reflines = setdiff(acslines, datlines);
  NRef = length(reflines);
  
  
  return;
  


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

  
  
  %toklines = zeros(1, NImg);
  %toklines(datlines) = toklines(datlines) + 1;
  %toklines(acslines) = toklines(acslines) + 1;
  %toklines(reflines) = toklines(reflines) + 1;
