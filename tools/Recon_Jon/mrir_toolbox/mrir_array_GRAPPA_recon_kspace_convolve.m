function [i_conv, k_conv] = mrir_array_GRAPPA_recon_kspace_convolve(k_data, G, evp, varargin)
%MRIR_ARRAY_GRAPPA_RECON_KSPACE_CONVOLVE
%
% [i_conv, k_conv] = mrir_array_GRAPPA_recon_kspace_convolve(raw, G, evp)

% TODO: add zero padding around boundary for kernel margins

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/15
% $Id: mrir_array_GRAPPA_recon_kspace_convolve.m,v 1.3 2009/02/15 23:31:31 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  t0 = clock;


  %==--------------------------------------------------------------------==%

  METHOD = 'conv2';

  if ( nargin >= 4 ),
    METHOD = lower(varargin{1});

    dstr = sprintf('computing reconstruction with method "%s"', METHOD);
    disp(sprintf('<i> [%s]: %s', mfilename, dstr));
  end;

  DO__MATRIX_STACKING = 0;


  %==--------------------------------------------------------------------==%

  NCol = mrir_ice_dimensions(k_data, 'col');
  NLin = mrir_ice_dimensions(k_data, 'lin');
  NCha = mrir_ice_dimensions(k_data, 'cha');
  NPar = mrir_ice_dimensions(k_data, 'par');

  flReadoutOSFactor = 2;

  if ( (evp.NLinMeas / evp.NAFLin) >= NLin ),
    error('input k-space data must include skipped lines');
  end;
  
  k_redu = mrir_image_pad_meas(k_data, evp);
  RawLin = mrir_ice_dimensions(k_data, 'lin');
  RawPar = mrir_ice_dimensions(k_data, 'par');


  %==--------------------------------------------------------------------==%

  t1 = clock;

  if ( isstruct(G) ),
    k = mrir_array_GRAPPA_conv_kernel(G);
  else,
    k = G;
  end;

  
  switch METHOD,
    %----------------------------------------------------------------------%
   case 'conv2',

    k_conv_cha = complex(zeros(NCol, evp.NImageLins, NCha, NCha));

    for trg = 1:NCha,
      for src = 1:NCha,
        k_conv_cha(:,:,src,trg) = conv2( k_redu(:,:,src), k(:,:,1,src,trg), 'same');
      end;
    end;

    k_conv = squeeze(sum(k_conv_cha, 3));


    %----------------------------------------------------------------------%
   case {'convmtx2', 'cconvmtx2'},

    k_redu = double(k_redu);

    [k_test, zpad_lin, zpad_par] = paditup(k_redu, G, evp);


    % note that these are different from G.Nsrcx and G.Nsrcy
    Nx = size(k, 1);
    Ny = size(k, 2);

    Nc = size(k_test, 1);
    Nr = size(k_test, 2);
    
    N_conv_full = [(Nc+Nx-1), (Nr+Ny-1)];


    if ( DO__MATRIX_STACKING ),

      for trg = 1:NCha,
	disp(sprintf('   stacking: target channel %02d...', trg));

	K = [];
	dat_vec = [];
	for src = 1:NCha,
	  K = cat(2, K, mrir_convmtx2( k(:,:,1,src,trg) , Nc, Nr) );
	end;

	k_conv_vec = K * reshape(k_redu, [], 1);
	k_conv_mat = reshape(k_conv_vec, N_conv_full);

%	k_conv(:,:,trg) = k_conv_mat(1+(Nx-1)/2:end-(Nx-1)/2, 1+(Ny-1)/2:end-(Ny-1)/2);
	k_conv(:,:,trg) = k_conv_mat;

      end;

    else,  %   DO__MATRIX_STACKING == 0

      for trg = 1:NCha,
	disp(sprintf('   processing: target channel %02d...', trg));

	for src = 1:NCha,
	  K = convmtx2( k(:,:,1,src,trg) , Nc, Nr);

	  k_conv_vec = K * reshape(k_test(:,:,src), [], 1);
	  k_conv_mat = reshape(k_conv_vec, N_conv_full);

	  k_conv_cha(:,:,src,trg) = k_conv_mat(1+(Nx-1)/2:end-(Nx-1)/2, 1+(Ny-1)/2+zpad_lin:end-(Ny-1)/2-zpad_lin(2));

	end;
      end;

      k_conv = squeeze(sum(k_conv_cha, 3));

    end;

   otherwise,
    error('invalid convolution method "%s"', METHOD);
  end;


  t2 = clock;
  runtime_seconds = etime(t2,t1);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total recon time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  %==--------------------------------------------------------------------==%


  i_conv = mrir_conventional_2d(k_conv);

  tN = clock;
  runtime_seconds = etime(tN,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total execution time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  return;




%**************************************************************************%
function [k_padresize, zpad_lin, zpad_par] = paditup(k_data, G, evp)

  NLinMeas  = evp.NLinMeas;
  NParMeas  = evp.NParMeas;
  NFirstLin = evp.NFirstLin;
  NFirstPar = evp.NFirstPar;
  RawLin    = evp.RawLin;
  RawPar    = evp.RawPar;


  %==--------------------------------------------------------------------==%
  %%% extract size parameters from input data

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
      ( k_framelines1(1)   + [min(min(kernels1(:, 1:Nkernel1)))-1]*R1 ) : ...
      R1 : ...
      ( k_framelines1(end) + [max(max(kernels1(:, 1:Nkernel1)))+0]*R1 ) ;

  k_extendlines2 = ...
      ( k_framelines2(1)   + [min(min(kernels2(:, 1:Nkernel2)))-1]*R2 ) : ...
      R2 : ...
      ( k_framelines2(end) + [max(max(kernels2(:, 1:Nkernel2)))+0]*R2 ) ;

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

  k_padresize = mrir_zeropad(k_data,      zpad_below, 'pre');
  k_padresize = mrir_zeropad(k_padresize, zpad_above, 'post');

  
  zpad_lin = [zpad_below(2), zpad_above(2)];
  zpad_par = [zpad_below(9), zpad_above(9)];
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_recon_kspace_convolve.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
