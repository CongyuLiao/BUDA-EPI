function varargout = mrir_noise_synthesize(dims, noisecov, varargin)
%MRIR_NOISE_SYNTHESIZE
%
% noise_color = mrir_noise_synthesize(dims, noisecov);


% see also: "mrir_array_stats_sim_data"

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/23
% $Id: mrir_noise_synthesize.m,v 1.2 2009/03/10 01:00:28 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncha = dims(mrir_DIM_CHA);

  if ( (nargin < 2) || isempty(noisecov) ),
    noisecov = complex(rand(Ncha, Ncha), rand(Ncha, Ncha));
  end;

  if ( Ncha ~= size(noisecov, 1) ),
    error('channel number mismatch');
  end;

  FLAG__VERIFY = 0;
  if ( nargin >= 3 ),
    FLAG__VERIFY = varargin{1};
  end;



  %==--------------------------------------------------------------------==%

  noise_white = complex(randn(dims), randn(dims));

  % calculate projection from cholesky decomposition
  P = chol(noisecov/2);

  noise_color = mrir_array_transform(noise_white, P');


  %==--------------------------------------------------------------------==%

  if ( FLAG__VERIFY ),

    noisecov_color = mrir_array_stats_matrix(noise_color, 'cov', 0);

    v0 = noisecov(      logical(triu(ones(size(noisecov)))));
    v1 = noisecov_color(logical(triu(ones(size(noisecov)))));

    v_mag_err_rel = ( (v0 - v1) ./ v1 );
    v_mag_err_rms = sqrt(mean(abs(v_mag_err_rel).^2));

    disp(sprintf('correlated data exhibits prescribed covariance with [[%2.1f%%]] accuracy', 100 * [1 - (1/v_mag_err_rms)]));

  end;

  
  %==--------------------------------------------------------------------==%

  if ( nargout >= 1 ),
    varargout{1} = noise_color;
  end;
  
  if ( nargout >= 2 ),
  
    noisecov_synth = mrir_array_stats_matrix(noise_color, 'cov', 0);
    varargout{2} = noisecov_synth;
    
  end;
  
    
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_noise_synthesize.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

