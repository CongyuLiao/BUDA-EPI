function noise = mrir_array_GRAPPA_noise_compute(noisecov, G, NLinMeas, NParMeas, varargin)
%

%% superseded by "mrir_array_GRAPPA_2d_variance"


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/07
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  
  % G:  [NCha * (R1-1) * (R2-1)] x [NCha * Nsrcx * Nsrcy * Nsrcz]
  
  %%% S_trg: [NCha * (R1-1) * (R2-1)] x Nfit
  %%% S_src: [NCha * Nsrcx * Nsrcy * Nsrcz] x Nfit
  
  
  Ntrg = size(G.kernel{1},1);
  NCha = Ntrg / max([1,G.R1-1]) / max([1,G.R2-1]);
  

  cha = 0;
  noise = [];
  for target = 1:Ntrg,
    
    if ( mod(target-1, (Ntrg)/NCha) == 0 ),
      cha = cha + 1;
      noise(end+1) = noisecov(cha,cha);
    end;
    
    % weights for all source points
    w = G.kernel{1}(target, :);
    
    % source weights reshaped so each column is contribution from one coil
    W = reshape(w, [], NCha).';

    % compute trace since only want to sum up for each source point (i.e., there is no correlation between different source points)
    noise(end+1) = trace(W'*noisecov*W);
    
  end;
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
