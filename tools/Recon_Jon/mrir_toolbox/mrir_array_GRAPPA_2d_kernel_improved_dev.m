function varargout = mrir_array_GRAPPA_2d_kernel_improved_dev(k_ACS, R1, R2, varargin)
%MRIR_ARRAY_GRAPPA_2D_KERNEL
%
% G = mrir_array_GRAPPA_2d_kernel(k_ACS, R1, R2)
% G = mrir_array_GRAPPA_2d_kernel(k_ACS, R1, R2, Nx, Ny, Nz)
% G = mrir_array_GRAPPA_2d_kernel(k_ACS, R1, R2, Nx, Ny, Nz, Nk1, Nk2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/10
% $Id: mrir_array_GRAPPA_2d_kernel.m,v 1.8 2009/10/03 21:47:15 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.8 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;
  global VERBOSE; if ( isempty(VERBOSE) ), VERBOSE = 1; end;

  global REGULARIZE; if ( isempty(REGULARIZE) ), REGULARIZE = 0; end;


  %==--------------------------------------------------------------------==%

  t0 = clock;


  Nsrcx = 3;        % kernel extent in frequency-encoded direction
  Nsrcy = 3;        % kernel extent in primary phase-encoded direction   (a.k.a., Nblocks1)
  Nsrcz = 3;        % kernel extent in secondary phase-encoded direction (a.k.a., Nblocks2)
  Nkernel1 = 1;     % number of shifted kernels: 1 <= Nkernel1 < Nblocks1
  Nkernel2 = 1;     % number of shifted kernels: 1 <= Nkernel2 < Nblocks2


  % kernel parameters

  if ( nargin >= 4 ),
    Nsrcx  = varargin{4-3};
  end;

  if ( nargin >= 5 ),
    Nsrcy  = varargin{5-3};
  end;

  if ( nargin >= 6 ),
    Nsrcz  = varargin{6-3};
  end;



  if ( nargin >= 7 ),
    Nkernel1  = varargin{7-3};
  end;

  if ( Nkernel1 > Nsrcy ),
    warning(sprintf('Nkernel1 set to maximum: %d', Nsrcy));
    Nkernel1 = Nsrcy;
  end;


  if ( nargin >= 8 ),
    Nkernel2  = varargin{8-3};
  end;

  if ( Nkernel2 > Nsrcz ),
    warning(sprintf('Nkernel2 set to maximum: %d', Nsrcz));
    Nkernel2 = Nsrcz;
  end;

  eta = 1; %normalization
  chi = 0.01; %  lower to deregularize which cause lower SNR
  %iExtRefColRange = NCol
  iExtRefColRange = 64;
  if ( nargin >= 9 ),
    eta = varargin{9-3}
    if ( eta > 5 ),
      warning(sprintf('value of eta (%f) seems too large (based on heuristic) -- consider using smaller value', eta));
    end;
  end;
 
  if ( nargin >= 10 ),
    chi = varargin{10-3}
  end;
  
  if ( nargin >= 11 ),
    iExtRefColRange = varargin{11-3}
  end;
  
  
  %==--------------------------------------------------------------------==%

  % the kernel sizes need to be included with the kernel itself since the
  % same values must be used during reconstruction. since they cannot be
  % deduced from the kernel size (i.e., multiple value pairs could yield the
  % same dimensions), for safety we adopt a convention where they are stored
  % in the kernel struct
  G.kernel = cell(Nkernel1, Nkernel2);
  G.Nsrcx = Nsrcx;
  G.Nsrcy = Nsrcy;
  G.Nsrcz = Nsrcz;

  G.R1 = R1;
  G.R2 = R2;


  NCol = mrir_ice_dimensions(k_ACS, 'col');
  NCha = mrir_ice_dimensions(k_ACS, 'cha');


  % number of collected ACS lines
  Nacs1 = mrir_ice_dimensions(k_ACS, 'lin');
  Nacs2 = mrir_ice_dimensions(k_ACS, 'par');

  if ( VERBOSE ),

%     mrir_array_GRAPPA_2d_kernel__display_params(NCol, NCha, ...
%                                                 Nacs1, Nacs2, ...
%                                                 Nsrcx, Nsrcy, Nsrcz, ...
%                                                 R1, R2, ...
%                                                 Nkernel1, Nkernel2);
% 

     mrir_array_GRAPPA_2d_kernel__display_params(min([NCol, 64]), NCha, ...
                                                 Nacs1, Nacs2, ...
                                                 Nsrcx, Nsrcy, Nsrcz, ...
                                                 R1, R2, ...
                                                 Nkernel1, Nkernel2);
  end;

  % TODO: add regularization and normalization parameters to display
  G.eta = eta;
  G.chi = chi;
  
  
  %==--------------------------------------------------------------------==%
  %%% generate all possible kernels given number of source points per kernel

  kernels1 = mrir_array_GRAPPA__pick_kernels(Nsrcy);
  kernels2 = mrir_array_GRAPPA__pick_kernels(Nsrcz);


  %==--------------------------------------------------------------------==%

  width_x = (Nsrcx-1)/2;
  kernel_x = ceil(-width_x : +width_x);


  iFirstValidCol = (NCol-iExtRefColRange)/2 + 1;
  iLastValidCol  = (NCol+iExtRefColRange)/2;
  
  k_interiorlines0 = iFirstValidCol:iLastValidCol;
  
  
  for ind1 = 1:Nkernel1,

    kernel_y = kernels1(:,ind1).' * R1;

    % for this kernel, only train on ACS lines
    k_interiorlines1 = 1 - min(kernel_y) : Nacs1 - max([(R1-1), max(kernel_y)]);

    for ind2 = 1:Nkernel2,

      kernel_z = kernels2(:,ind2).' * R2;
      k_interiorlines2 = 1 - min(kernel_z) : Nacs2 - max([(R2-1), max(kernel_z)]);


      % number of fits of the GRAPPA kernel obtained from this data
      Nfit = length(k_interiorlines0) * length(k_interiorlines1) * length(k_interiorlines2);

%      fprintf(1, '<i> [%s]: allocating...', mfilename);
      % preallocate for time efficiency
      S_src = zeros(NCha * Nsrcx * Nsrcy * Nsrcz, Nfit);
%      fprintf(1, 'and...');
      S_trg = zeros(NCha * (R1*R2-1),             Nfit);
%      fprintf(1, 'done!\n');

      t1 = clock;

      for kk = 1:length(k_interiorlines2),
        for jj = 1:length(k_interiorlines1),
          for ii = 1:length(k_interiorlines0),

            % index of current loop iteration
            fitnum = sub2ind([length(k_interiorlines0), length(k_interiorlines1), length(k_interiorlines2)], ii, jj, kk);

            % "pos" designates center of kernel
            Xpos = k_interiorlines0(ii);

            Ypos = k_interiorlines1(jj);
            Zpos = k_interiorlines2(kk);


            Xsrc = mod([ (Xpos-1) + kernel_x ], NCol) + 1;

            % to enable wrap-around:
            %    Ysrc = mod([ (Ypos-1) + kernel_y ], Nacs1) + 1;
            %    Zsrc = mod([ (Zpos-1) + kernel_z ], Nacs2) + 1;

            Ysrc = Ypos + kernel_y;
            Zsrc = Zpos + kernel_z;


            % ordering convention is that coil channels are at outside of loop
            s_src = mrir_array_GRAPPA__extract_entries(k_ACS, Xsrc, Ysrc, Zsrc, 1:NCha);


            [Xtrg, Ytrg, Ztrg] = ndgrid( Xpos, Ypos:(Ypos+(R1-1)), Zpos:(Zpos+(R2-1)) );
            if ( length(Xtrg) > 1 ),
              Xtrg(1) = [];
              Ytrg(1) = [];
              Ztrg(1) = [];
            end;

            %    Ytrg = mod([ (Ytrg-1) ], Nacs1) + 1;
            %    Ztrg = mod([ (Ztrg-1) ], Nacs2) + 1;


            % force slowest looping over channels (as with source point indices)
            [xtrg, ctrg] = ndgrid(Xtrg, 1:NCha);
            [ytrg, ctrg] = ndgrid(Ytrg, 1:NCha);
            [ztrg, ctrg] = ndgrid(Ztrg, 1:NCha);


            i_trg = sub2ind(size(squeeze(k_ACS)), xtrg(:), ytrg(:), ctrg(:), ztrg(:));

            s_trg = k_ACS(i_trg);


            %%% S_src: [NCha * Nsrcx * Nsrcy * Nsrcz] x Nfit
            S_src(:, fitnum) = s_src;


            %%% S_trg: [NCha * (R1 * R2 - 1)] x Nfit
            S_trg(:, fitnum) = s_trg;

          end;  %% ii
        end;  %% jj
      end;  %% kk

      t2 = clock;

      if ( VERBOSE ),
      setuptime_seconds = etime(t2,t1);
      TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                     fix(setuptime_seconds/60/60), ...
                     rem(fix(setuptime_seconds/60), 60), ...
                     rem(fix(setuptime_seconds), 60));

      dstr = sprintf('    setup time = %s', TIME);
      disp(sprintf('<t> [%s]: %s', mfilename, dstr));


      disp(sprintf('==> [%s]: density of source matrix (should be near 1.00) = %2.3f', mfilename, nnz(S_src) / numel(S_src)));


      dstr = sprintf('computing pseudo-inverse of %dx%d source matrix...', size(S_src,1), size(S_src,2));
      disp(sprintf('<i> [%s]: %s', mfilename, dstr));
      end;


      %==--------------------------------------------------------------------==%
      % G:  [NCha * (R1 * R2 - 1)] x [NCha * Nsrcx * Nsrcy * Nsrcz]

      % note that if either eta == 0 or chi == 0 naturally their respective
      % normalization/regularization is disabled
      
      % "ArtifactReduction": default eta = 1.0
	
	colskip = length(k_interiorlines0);
	
	Psi = [];
	for col = 1:colskip:fitnum,

	  bcol = S_trg(:,col:(col+colskip-1));

	  p = sum(sum(abs(bcol).^2));
	  Psi(end+1) = p.^(-eta/2);

	  if ( Psi < sqrt(eps) ),
	    error('value of eta (%f) seems too large -- numerical truncation detected', eta);
	  end;
	  
	  S_trg(:,col:(col+colskip-1)) = S_trg(:,col:(col+colskip-1)) * Psi(end);
	  S_src(:,col:(col+colskip-1)) = S_src(:,col:(col+colskip-1)) * Psi(end);

	end;

      
      % "NoiseReductionI": default chi = 0.0001

      % Ghat = ( B A' ) ( A' A + lambda E)^-1 
      
      U = S_trg * S_src';

      M = S_src * S_src';
      t = trace(M);

      lambda = chi * t / size(S_src, 1);

      M = M + ( lambda * eye(size(M)) );

      V = inv(M);


      % store kernel computed for this block
      G.kernel{ind1, ind2} = U * V;


      condition_number = cond(G.kernel{ind1, ind2});


      %==--------------------------------------------------------------------==%

      if ( VERBOSE ),
        dstr = sprintf('condition number of %dx%d source matrix:   ( %7.4f )', size(S_src,1), size(S_src,2), condition_number);
        disp(sprintf('<i> [%s]: %s', mfilename, dstr));
      end;

      % store condition number in kernel struct
      G.cond{ind1, ind2} = condition_number;

      % store vector norm in kernel struct
      G.norm{ind1, ind2} = norm( G.kernel{ind1, ind2}(:) );

      % allocate field for later use
      G.MSE{ind1, ind2} = [];

      t3 = clock;

      if ( VERBOSE ),
      pinvtime_seconds = etime(t3,t2);
      TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                     fix(pinvtime_seconds/60/60), ...
                     rem(fix(pinvtime_seconds/60), 60), ...
                     rem(fix(pinvtime_seconds), 60));

      dstr = sprintf('   "pinv" time = %s', TIME);
      disp(sprintf('<t> [%s]: %s', mfilename, dstr));
      end;
    end;
  end;


  %==--------------------------------------------------------------------==%

  tN = clock;
  runtime_seconds = etime(tN,t0);

  G.computed = datestr(tN, 'yyyy-mmm-dd HH:MM:SS');

  if ( VERBOSE ),
  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total run time = %s', TIME);
  disp(sprintf('\n<t> [%s]: %s', mfilename, dstr));
  end;

  if ( nargout > 0 ),
    varargout{1} = G;
  end;


  return;



%**************************************************************************%
function [m, X, Y, Z, C] = mrir_array_GRAPPA__extract_entries(M, x, y, z, c)

% extract entries from volume based on matrices of coordinates

  [X, Y, C, Z] = ndgrid(x, y, c, z);
  ind = sub2ind(size(M), X, Y, C, Z);
  m = reshape(M(ind), [], 1);


  return;



%**************************************************************************%
function mrir_array_GRAPPA_2d_kernel__display_params(NCol, NCha, Nacs1, Nacs2, ...
                                                  Nsrcx, Nsrcy, Nsrcz, R1, R2, ...
                                                  Nkernel1, Nkernel2)

  dstr = sprintf('<i> [%s]: GRAPPA kernel fitting parameters...\n', mfilename);

  dstr = strvcat(dstr, sprintf('\t       ACS dimensions:  NCol = %3d, NLin = %3d, NPar = %3d, NCha = %3d', NCol, Nacs1, Nacs2, NCha));
  dstr = strvcat(dstr, sprintf('\t    kernel dimensions:  Nsrcx = %2d, Nsrcy = %2d, Nsrcz = %2d', Nsrcx, Nsrcy, Nsrcz));
  dstr = strvcat(dstr, sprintf('\t                        number of target points per channel = %d', (R1*R2)-1));
  dstr = strvcat(dstr, sprintf('\t block reconstruction:  Nkernel1 =%2d, Nkernel2 =%2d', Nkernel1, Nkernel2));


  disp(dstr);
  fprintf(1, '\n');

  return;


%**************************************************************************%
function kernels = mrir_array_GRAPPA__pick_kernels(Nsrc)


  kernels = [ones(Nsrc,1)  *     (0:Nsrc-1)] ...
            - [(0:Nsrc-1)' * ones(1,Nsrc)];

  % a lazy way to re-order the kernels inside-out:
  [sorted, col_permute] = sort(abs( ([-(Nsrc-1)/2]:[+(Nsrc-1)/2]) - 0.25 ) );
  kernels = sortrows(kernels(:,col_permute));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_2d_kernel.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 85
  %%% comment-column: 0
  %%% End:
