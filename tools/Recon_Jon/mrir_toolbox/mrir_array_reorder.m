function [data_sort, channel_ind, elements_sorted]  = mrir_array_reorder(data, varargin)
%MRIR_ARRAY_REORDER
%
% data_sort = mrir_array_reorder(data, element_str)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/07
% $Id: mrir_array_reorder.m,v 1.2 2008/11/07 08:02:02 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( iscell(varargin{1}) ),
    element_str = varargin{1};
    [elements_sorted, channel_ind] = sort(element_str);
  else,
    channel_ind = varargin{1};
  end;
  
  
  if (  ( ndims(data) == 2 ) && ( length(element_str) == size(data,1) ) && ( length(element_str) == size(data,2) )  ),

    data_sort    = data(channel_ind,channel_ind);
  
  else,
    %                   1 2            3  4 5 6 7 8 9 0 1 2 3 4 5 6            
    data_sort    = data(:,:, channel_ind, :,:,:,:,:,:,:,:,:,:,:,:,:);
  end;
  


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_reorder.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


