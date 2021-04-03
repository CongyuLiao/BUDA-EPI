function [img_map_flip, img_synth_volt, img_synth_flip, img_map_volt, cost, param1, param2] = mrir_transmit_map(img, voltages, varargin)
%MRIR_TRANSMIT_MAP
%
% img_map_flip = mrir_transmit_map(img, voltages, desired_volt)
%
% [img_map_flip, img_synth_volt, img_synth_flip, img_map_volt] = ...
%      mrir_transmit_map(img, voltages, desired_volt, desired_flip)
%
%
% images acquired across a range of transmit voltages should be concatenated
% along dimension 7 ("mrir_DIM_REP") of the image data array.
%
% example:
%
%  img_FA = cat(7, ...
%               img_FA_voltage1, ...
%               img_FA_voltage2, ...
%               img_FA_voltage3, ...
%               img_FA_voltage4, ...
%               img_FA_voltage5);
%
%  img_map_flip = mrir_transmit_map(img_FA, voltages, prot.flAmplitude, prot.adFlipAngleDegree/180*pi);
%
%
% all input and output flip angles are assumed to be in radians.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/12
% $Id: mrir_transmit_map.m,v 1.4 2008/10/23 16:30:47 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%
  %%% set defaults to be typical parameter values

  desired_volt = 150;
  desired_flip = 90    / 180 * pi;

  if ( nargin >= 3 ),
    desired_volt = varargin{1};
  end;

  if ( nargin >= 4 ),
    desired_flip = varargin{2};
  end;


  %==--------------------------------------------------------------------==%

  voltages_data = voltages(:);
  Nvoltages = length(voltages);

  [voltages_data, sort_ind] = sort(voltages_data);

  dims = size(img);
  dims(end+1:16) = 1;

  %                1 2 3 4 5 6        7 8 9 0 1 2 3 4 5 6
  img_matrix = img(:,:,:,:,:,:,sort_ind,:,:,:,:,:,:,:,:,:);

  % rearrange image data so acquisitions across voltages are first
  perm = [7, 1:6, 8:16];
  img_vector = reshape(permute(img, perm), dims(7), []);

  % options for optimization routine
  options = optimset(@fminsearch);
  options.Display     = 'off';
  options.MaxFunEvals = 50000;
  options.MaxIter     = 10000;

  img_max = max(img_vector(:));

  % make an educated guess for the parameter values to start off the fitting
  X0 = [img_max/2; 1/(50*2*pi)];

  % for fit of first pixel, initialize optimization with the guess, and for
  % subsequent pixels use previously computed optimal parameter values
  % (under the assumption that the values vary smoothly, while the
  % optimization will still converge even if initialization is a little off)
  Xinit = X0;

  % allocate
  param    = zeros(2, size(img_vector,2));
  cost     = zeros(1, size(img_vector,2));

  if ( DEBUG ),
    fig = figure(100); figure(100); set(fig, 'Name', mfilename); ax = axis; box on;
  end;

  for ind = 1:size(img_vector,2),

    intensities = double(img_vector(:,ind));

    % skip if this pixel is all zeros (i.e., its been masked out)
    if ( ~any(intensities) ),
      continue;
    end;

    % HACK #1: heuristic to include enough points for accurate fit
    flip = find(diff(intensities(3:end) > intensities(1:end-2)) == -1);
    if ( isempty(flip) ),
      N = Nvoltages;
    else,
      N = flip(1) + 3;
    end;

    [X,fval,exitflag,output] = fminsearch(@mrir_transmit_map__flipangle_model_fit, ...
                                          Xinit, options, voltages_data(1:N), intensities(1:N));

    param(:,ind) = X(:);
    cost(:,ind) = fval;

    % HACK #2: only keep result to initialize next pixel if error is less than an arbitrary value
    if ( fval < 0.20 ),
      Xinit = X;
    end;


    if ( DEBUG ),
      voltages_plot = union(linspace(0, max(voltages_data), 20), voltages_data).';
      cla; hold on;
      set(gca, 'XLim', [0, round(voltages_plot(end)*1.1)], 'YLim', [0, round(img_max*1.1)]);
      plot(voltages_plot, mrir_transmit_map__flipangle_model(X, voltages_plot), 'r-');
      plot(voltages_data(1:N), intensities(1:N), 'bo');
      plot(desired_volt,  mrir_transmit_map__flipangle_model(X, desired_volt),  'gx')
      [r, c] = ind2sub([size(img, 1), size(img, 2)], ind);
      title(sprintf('(%03d,%03d) cost: %2.2f, N=%2d', r, c, cost(:,ind)*100, N));
      pause(0.001);
    end;
  end;

  % synthetic image for a given voltage V0:            I = A * sin ( V0 * B )
  img_synth_volt = mrir_transmit_map__flipangle_model(param, desired_volt);

  % synthetic image for a given flip angle theta0:     I = A * sin ( theta0 )
  img_synth_flip = param(1,:) .* sin(desired_flip);

  % mapping of flip angle for a given voltage V0:      M = V0 * B
  img_map_flip = desired_volt .* param(2,:);

  % mapping of voltage for a given flip angle theta0:  M = theta0 / B
  img_map_volt = desired_flip ./ param(2,:);


  dims = [1, dims(1:6), dims(8:16)];

  img_synth_volt = ipermute(reshape(img_synth_volt, dims), perm);
  img_synth_flip = ipermute(reshape(img_synth_flip, dims), perm);

  img_map_flip = ipermute(reshape(img_map_flip, dims), perm);
  img_map_volt = ipermute(reshape(img_map_volt, dims), perm);


  param1 = ipermute(reshape(param(1,:), dims), perm);
  param2 = ipermute(reshape(param(2,:), dims), perm);
  cost   = ipermute(reshape(cost, dims), perm);


  return;



%**************************************************************************%
function estimate = mrir_transmit_map__flipangle_model(param, voltages)

%  y = A * cos( B*x )
%      y: intensity
%      x: voltages
%      A: amplitude
%      B: period


  A = param(1, :);
  B = param(2, :);

  estimate = A .* sin( voltages * B );


  return;



%**************************************************************************%
function cost = mrir_transmit_map__flipangle_model_fit(param, voltages, intensities)

  estimate = mrir_transmit_map__flipangle_model(param, voltages);

  residual = (estimate - intensities) ./ mean(intensities);

  % RMS
  cost = sqrt(mean(residual.^2));

  % MED
  %%%  cost = median(residual);

  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_transmit_map.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
