function varargout = mrir_array_GRAPPA_2d_variance(noisecov, k_reduFOV, G, varargin)
%MRIR_ARRAY_GRAPPA_2D_VARIANCE
%
% v_fullFOV = mrir_array_GRAPPA_2d_variance(noisecov, k_reduFOV, G, NLinMeas, NParMeas)
%
% v_fullFOV = mrir_array_GRAPPA_2d_variance(noisecov, k_reduFOV, G, NLinMeas, NParMeas, NFirstLin, NFirstPar)
%
% v_fullFOV = mrir_array_GRAPPA_2d_variance(noisecov, k_reduFOV, G, evp)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/10
% $Id: mrir_array_GRAPPA_2d_variance.m,v 1.1 2008/12/09 22:13:54 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  t0 = clock;
  
  % data collection parameters:
  
  
  if ( isstruct(varargin{1}) ),
    evp = varargin{1};
    NLinMeas = evp.NLinMeas;
    NParMeas = evp.NParMeas;

    NFirstLin = evp.NFirstLin;
    NFirstPar = evp.NFirstPar;
  
  else,
  
    NLinMeas = varargin{1};
    NParMeas = varargin{2};
    
    NFirstLin = 1;
    NFirstPar = 1;

    if ( nargin >= 6 ),
      NFirstLin = varargin{3};
    end;
    
    if ( nargin >= 7 ),
      NFirstPar = varargin{4};
    end;

  end;  


  %==--------------------------------------------------------------------==%
  %%% extract size parameters from input data

  dims_redu = size(k_reduFOV);

  NCol = mrir_ice_dimensions(k_reduFOV, 'col');
  NCha = mrir_ice_dimensions(k_reduFOV, 'cha');
  RawLin = mrir_ice_dimensions(k_reduFOV, 'lin');
  RawPar = mrir_ice_dimensions(k_reduFOV, 'par');

  [Nkernel1, Nkernel2] = size(G.kernel);

  % kernel parameters:
  Nsrcx = G.Nsrcx;  % kernel extent in frequency-encoded direction       (odd)
  Nsrcy = G.Nsrcy;  % kernel extent in primary phase-encoded direction   (even)
  Nsrcz = G.Nsrcz;  % kernel extent in secondary phase-encoded direction (even)

  R1 = G.R1;  % acceleration factor in frequency-encoded direction
  R2 = G.R2;  % acceleration factor in primary phase-encoded direction


  %==--------------------------------------------------------------------==%
  %%% generate all possible kernels given number of source points per kernel

  kernels1 = mrir_array_GRAPPA__pick_kernels(Nsrcy);
  kernels2 = mrir_array_GRAPPA__pick_kernels(Nsrcz);


  %==--------------------------------------------------------------------==%
  %%% calculate lines in k-space containing data (k-space-based indexing)

  % its AMAZING that this works...
  k_datlines1 = NFirstLin : R1 : NLinMeas;
  full2dat1(k_datlines1) = 1:length(k_datlines1);
  if ( length(k_datlines1) ~= RawLin ), warning('i am confused! (1)'); end;

  k_datlines2 = NFirstPar : R2 : NParMeas;
  full2dat2(k_datlines2) = 1:length(k_datlines2);
  if ( length(k_datlines2) ~= RawPar ), warning('i am confused! (2)'); end;

  % "NLastLin" and "NLastPar" are the indices of last lines that actually
  % contain data
  NLastLin = k_datlines1(end);  % previously: RawLin * R1
  NLastPar = k_datlines2(end);  % previously: RawPar * R2


  %  % for some siemens data, RawLin*R1 <= NLinMeas, but by an amount less than
  %  % R.
  %  NLinMeas == (RawLin*R1 + NFirstLin-1)
  %  NParMeas == (RawPar*R2 + NFirstPar-1)
  %
  %  if ( datmax1 ~= (RawLin * R1) ), warning('i am confused! (3)'); end;
  %  if ( datmax2 ~= (RawPar * R2) ), warning('i am confused! (4)'); end;
  %
  %  %(see, e.g., "meas_MID714_GRE3D_grappa_2x2_int_FID9798")


  %==--------------------------------------------------------------------==%
  %%% TODO: evaluate cyclic indices to handle boundaries witout zero padding

  % wrap condition 1: (mod( [datmax1 + R1 - 1], NLinMeas ) + 1) == NFirstLin
  % wrap condition 2: (mod( [datmax2 + R2 - 1], NParMeas ) + 1) == NFirstPar


  %==--------------------------------------------------------------------==%
  %%% augment k-space array with zeros before applying kernel


  % calculate lines needed to frame k-space so that we can zero pad data
  % lines to surround all missing lines.
  k_framelines1 = [1 + rem(NFirstLin - R1 - 1, R1)] : R1 : [NLinMeas + rem(NLastLin + R1 - NLinMeas, R1)];
  k_framelines2 = [1 + rem(NFirstPar - R2 - 1, R2)] : R2 : [NParMeas + rem(NLastPar + R2 - NParMeas, R2)];

  % calculate how many lines are needed in order to further surround data so
  % that whenever the kernel is centered on a missing data line it extends
  % into a margin of zeros. (here we need the max and min across all kernels
  % since we will be averaging together reconstructions across kernels.)
  k_extendlines1 = ...
      ( k_framelines1(1)   + [min(min(kernels1(:, 1:Nkernel1)))+0]*R1 ) : ...
      R1 : ...
      ( k_framelines1(end) + [max(max(kernels1(:, 1:Nkernel1)))-1]*R1 ) ;

  k_extendlines2 = ...
      ( k_framelines2(1)   + [min(min(kernels2(:, 1:Nkernel2)))+0]*R2 ) : ...
      R2 : ...
      ( k_framelines2(end) + [max(max(kernels2(:, 1:Nkernel2)))-1]*R2 ) ;

  if ( Nsrcz == 1 ),
    k_extendlines2 = 1;
  end;

  k_ind_shift1 = 1 - k_extendlines1(1);
  k_ind_shift2 = 1 - k_extendlines2(1);

  % given the extents, calculate how much we need to zero pad data. (data has 16 dimensions)
  [zpad_below, zpad_above] = deal(zeros(16,1));

  % padding needed for LIN dimension
  zpad_below(2) = k_ind_shift1;
  zpad_above(2) = k_extendlines1(end) - NLastLin;

  % padding needed for PAR dimension
  zpad_below(9) = k_ind_shift2;
  zpad_above(9) = k_extendlines2(end) - NLastPar;

  if ( Nsrcz == 1 ),
    zpad_above(9) = 0;
  end;

  % fill empty array with data lines at the appropriate k-space indices
  % (pre-allocation needed for large data sets)
  dims_full = ones(1,16);
  dims_full(1:length(dims_redu)) = dims_redu;
  dims_full(2) = NLinMeas;
  dims_full(9) = NParMeas;

  k_padresize = zeros(dims_full, 'single');
  %           1           2 3 4 5 6 7 8           9
  k_padresize(:,k_datlines1,:,:,:,:,:,:,k_datlines2) = repmat(reshape(diag(noisecov), 1, 1, NCha), [NCol, length(k_datlines1), 1, 1,1,1,1,1, length(k_datlines2)]);

  k_padresize = padarray(k_padresize, zpad_below, 'pre');
  k_padresize = padarray(k_padresize, zpad_above, 'post');


  % convert k-space positions into (1-based) array indices by index shift
  i_framelines1 = k_framelines1 + k_ind_shift1;
  i_framelines2 = k_framelines2 + k_ind_shift2;

  % HACK for 1D reconstructions
  if ( Nsrcz == 1 ),
    i_framelines2 = [1 1];
  end;

  %==--------------------------------------------------------------------==%

  %%% TODO: account for positions on end of k-space where padding results in
  %%% 0-variance data points to be mixed in, reducing the variance
  
  
  Ntrg = size(G.kernel{1},1);

  for target = 1:Ntrg,
    
    % weights for all source points
    w = G.kernel{1,1}(target, :);
    
    % source weights reshaped so each column is contribution from one coil
    W = reshape(w, [], NCha).';

    % compute trace since only want to sum up across coils for each source
    % point (i.e., there is no correlation between non-corresponding source
    % points)
    noise(target) = trace(W'*noisecov*W);
    
  end;
  

  %==--------------------------------------------------------------------==%

      for Xkern = 1:NCol,
        for Ykern = i_framelines1(1:end-1),
          for Zkern = i_framelines2(1:end-1),

            % source and target positions all in terms of (1-based) array
            % indices!

            % grab target points above current data line, between current and next
            [Xtrg, Ytrg, Ztrg] = ndgrid( Xkern, Ykern:(Ykern+(R1-1)), Zkern:(Zkern+(R2-1)) );

            % remove first point only, which will correspond to a source
            % point; this step also converts arrays into vectors.
            Xtrg(1) = [];
            Ytrg(1) = [];
            Ztrg(1) = [];

            %    Ytrg = mod([ (Ytrg-1) ], NLinMeas) + 1;
            %    Ztrg = mod([ (Ztrg-1) ], NParMeas) + 1;

	    % force slowest looping over channels (as with source point indices)
	    [xtrg, ctrg] = ndgrid(Xtrg, 1:NCha);
            [ytrg, ctrg] = ndgrid(Ytrg, 1:NCha);
            [ztrg, ctrg] = ndgrid(Ztrg, 1:NCha);

	    
            i_trg = sub2ind(size(squeeze(k_padresize)), xtrg(:), ytrg(:), ctrg(:), ztrg(:));
	    
            % insert target data points into array; this step takes MOST of
            % the computation time!
            k_padresize(i_trg) = noise;

          end;
        end;
      end;


  %==--------------------------------------------------------------------==%
  %%% extract out only data lines requested by caller from augmented/padded array

  %                       1                         2 3 4 5 6 7 8                         9 0 1 2 3 4 5 6
  k_fullFOV = k_padresize(:,k_ind_shift1+(1:NLinMeas),:,:,:,:,:,:,k_ind_shift2+(1:NParMeas),:,:,:,:,:,:,:);


  if ( nargout > 0 ),
    varargout{1} = k_fullFOV;
  end;


  %==--------------------------------------------------------------------==%

  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total computation time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  return;



%**************************************************************************%
function [m, X, Y, Z, C, ind] = mrir_array_GRAPPA__extract_entries(M, x, y, z, c)



  [X, Y, C, Z] = ndgrid(x, y, c, z);
  ind = sub2ind(size(M), X, Y, C, Z);
  ind = permute(ind, [1, 2, 4, 3]);
 
  m = M(ind(:));


  return;



%**************************************************************************%
function M = mrir_array_GRAPPA__insert_entries(M, v, x, y, z, c, N1, N2)


  [X, Y, C, Z] = ndgrid(x, y, c, z);
  ind = sub2ind(size(squeeze(M)), X, Y, C, Z);
  ind = permute(ind, [1, 2, 4, 3]);
  
  M(ind) = M(ind(:)) + (v/N1/N2);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_variance.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
