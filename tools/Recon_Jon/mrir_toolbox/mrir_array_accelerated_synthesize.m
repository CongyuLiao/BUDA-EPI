function varargout = mrir_array_accelerated_synthesize(meas, evp, FLAG__KSPACE)
%MRIR_ARRAY_ACCELERATED_SYNTHESIZE
%
% mrir_array_accelerated_synthesize(dat, evp, FLAG__KSPACE)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/23
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  datlines = evp.NFirstLin : evp.NAFLin : evp.NLinMeas;

  if ( length(datlines) ~= evp.RawLin ),
    error('number of extracted lines does not match number in header');
  end;


  acslines = evp.NFirstRefLin : (evp.NFirstRefLin + evp.NRefLin - 1);


  %          1        2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  dat = meas(:,datlines,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

  acs = meas(:,acslines,:,:,:,:,:,:,:,:,:,:,:,:,:,:);



  if ( nargout > 0 ),

    varargout{1} = dat;
    varargout{2} = acs;

  end;


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
