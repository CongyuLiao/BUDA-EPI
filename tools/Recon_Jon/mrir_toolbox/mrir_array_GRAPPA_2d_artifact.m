function varargout = mrir_array_GRAPPA_2d_artifact(k_fullFOV, G, R1, R2, NLinMeas, NParMeas, varargin)
%MRIR_ARRAY_GRAPPA_2D_ARTIFACT
%
% k_artifact = mrir_array_GRAPPA_2d_artifact(k_fullFOV, G, R1, R2, NLinMeas, NParMeas)
%
% k_artifact = mrir_array_GRAPPA_2d_artifact(k_fullFOV, G, R1, R2, NLinMeas, NParMeas, NFirstLin, NFirstPar)


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/oct/19
% $Id: mrir_array_GRAPPA_2d_artifact.m,v 1.1 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  t0 = clock;


  NFirstLin = 1;
  NFirstPar = 1;


  % data collection parameters:

  if ( nargin >= 7 ),
    NFirstLin = varargin{7-6};
  end;

  if ( nargin >= 8 ),
    NFirstPar = varargin{8-6};
  end;


  %==--------------------------------------------------------------------==%
  %%% extract size parameters from input data

  dims_redu = size(k_fullFOV);
  
  NCol = mrir_ice_dimensions(k_fullFOV, 'col');
  NCha = mrir_ice_dimensions(k_fullFOV, 'cha');
  RawLin = mrir_ice_dimensions(k_fullFOV, 'lin');
  RawPar = mrir_ice_dimensions(k_fullFOV, 'par');

  [Nkernel1, Nkernel2] = size(G.kernel);

  % kernel parameters:
  Nsrcx = G.Nsrcx;  % kernel extent in frequency-encoded direction       (odd)
  Nsrcy = G.Nsrcy;  % kernel extent in primary phase-encoded direction   (even)
  Nsrcz = G.Nsrcz;  % kernel extent in secondary phase-encoded direction (even)


  %==--------------------------------------------------------------------==%
  %%% generate all possible kernels given number of source points per kernel
  
  kernels1 = mrir_array_GRAPPA__pick_kernels(Nsrcy);
  kernels2 = mrir_array_GRAPPA__pick_kernels(Nsrcz);
  

  %==--------------------------------------------------------------------==%
  %%% calculate lines in k-space containing data (k-space-based indexing)

  % its AMAZING that this works...
  k_datlines1 = NFirstLin : R1 : NLinMeas;
  k_datlines2 = NFirstPar : R2 : NParMeas;

  % "NLastLin" and "NLastPar" are the indices of last lines that actually
  % contain data
  NLastLin = k_datlines1(end);  % previously: RawLin * R1
  NLastPar = k_datlines2(end);  % previously: RawPar * R2


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

  k_augmented = k_fullFOV;

  k_augmented = padarray(k_augmented, zpad_below, 'pre');
  k_augmented = padarray(k_augmented, zpad_above, 'post');

  
  % convert k-space positions into (1-based) array indices by index shift
  ind_framelines1 = k_framelines1 + k_ind_shift1;
  ind_framelines2 = k_framelines2 + k_ind_shift2;
  
  if ( Nsrcz == 1 ),
    ind_framelines2 = [1 1];
  end;
  
  %==--------------------------------------------------------------------==%

  % initialize array for storing recons for each individual kernel
%  k_accumulate = zeros(size(k_augmented));
  
  width_x = (Nsrcx-1)/2;
  kernel_x = ceil(-width_x : +width_x);

  for ind1 = 1:Nkernel1,

    kernel_y = kernels1(:,ind1).' * R1;

    for ind2 = 1:Nkernel2,

      kernel_z = kernels2(:,ind2).' * R2;

      % copy the augmented k-space into new buffer for this kernel reconstruction
      k_buffer = k_augmented;
      
      
      for Xkern = 1:NCol,
        for Ykern = ind_framelines1(1):ind_framelines1(end-1),
          for Zkern = ind_framelines2(1):ind_framelines2(end-1),

	    % source and target positions all in terms of (1-based) array
            % indices!

            Xsrc = mod([ (Xkern-1) + kernel_x ], NCol) + 1;

            Ysrc = ( Ykern + kernel_y );
            Zsrc = ( Zkern + kernel_z );

            %    Ysrc = mod([ (Ykern-1) + kernel_y ], NLin) + 1;
            %    Zsrc = mod([ (Zkern-1) + kernel_z ], NPar) + 1;


            % grab target points above current data line, between current and next
            [Xtrg, Ytrg, Ztrg] = ndgrid( Xkern, Ykern:(Ykern+(R1-1)), Zkern:(Zkern+(R2-1)) );

            % remove first point only, which will correspond to a source
            % point; this step also converts arrays into vectors.
	    Xtrg(1) = [];
            Ytrg(1) = [];
            Ztrg(1) = [];

            %    Ytrg = mod([ (Ytrg-1) ], NLinMeas) + 1;
            %    Ztrg = mod([ (Ztrg-1) ], NParMeas) + 1;

            [xtrg, ctrg] = ndgrid(Xtrg, 1:NCha);
            [ytrg, ctrg] = ndgrid(Ytrg, 1:NCha);
            [ztrg, ctrg] = ndgrid(Ztrg, 1:NCha);

            i_trg = sub2ind(size(squeeze(k_augmented)), xtrg(:), ytrg(:), ctrg(:), ztrg(:));

            s_src = mrir_array_GRAPPA__extract_entries(k_augmented, Xsrc, Ysrc, Zsrc, 1:NCha);

            % apply kernel to source points to compute target points
	    s_trg = G.kernel{ind1, ind2} * s_src;

            % insert target data points into array; this step takes MOST of
            % the computation time!
            k_buffer(i_trg) = s_trg;

          end;
        end;
      end;

      % average in the results from reconstructing with this kernel
%      k_accumulate = k_accumulate + k_buffer;
      %k_buffer = [];
      
    end;
  end;

  %k_accumulate = k_accumulate / Nkernel1 / Nkernel2;
  

  %==--------------------------------------------------------------------==%
  %%% extract out only data lines requested by caller from augmented array

  %                        1                         2 3 4 5 6 7 8                          9 0 1 2 3 4 5 6
%  k_fullFOV = k_accumulate(:,k_ind_shift1+(1:NLinMeas),:,:,:,:,:,:,k_ind_shift2+(1:NParMeas),:,:,:,:,:,:,:);
  k_validate = k_buffer(:,k_ind_shift1+(1:NLinMeas),:,:,:,:,:,:,k_ind_shift2+(1:NParMeas),:,:,:,:,:,:,:);
%  k_accumulate = k_accumulate(:,k_ind_shift1+(1:NLinMeas),:,:,:,:,:,:,k_ind_shift2+(1:NParMeas),:,:,:,:,:,:,:);

  
  
  if ( nargout > 0 ),
    varargout{1} = k_validate;
%    varargout{2} = k_accumulate;
  end;


  %==--------------------------------------------------------------------==%

  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total recon time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  return;



%**************************************************************************%
function [m, X, Y, Z, C] = mrir_array_GRAPPA__extract_entries(M, x, y, z, c)


  [X, Y, C, Z] = ndgrid(x, y, c, z);
  ind = sub2ind(size(M), X, Y, C, Z);
  m = reshape(M(ind), [], 1);


  return;



%**************************************************************************%
function M = mrir_array_GRAPPA__insert_entries(M, v, x, y, z, c, N1, N2)


  [X, Y, C, Z] = ndgrid(x, y, c, z);
  ind = sub2ind(size(squeeze(M)), X, Y, C, Z);
  M(ind) = M(ind(:)) + (v/N1/N2);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_artifact.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
