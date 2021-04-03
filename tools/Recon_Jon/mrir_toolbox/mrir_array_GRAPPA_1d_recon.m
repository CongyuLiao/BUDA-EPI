function varargout = mrir_array_GRAPPA_1d_recon(k_reduFOV, G, R, Nlin, varargin)
%MRIR_ARRAY_GRAPPA_1D_RECON
%
% k_fullFOV = mrir_array_GRAPPA_1d(k_reduFOV, G, R, Nlin)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/10
% $Id: mrir_array_GRAPPA_1d_recon.m,v 1.2 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;
  
 
  %==--------------------------------------------------------------------==%

  t0 = clock;

  Nsrc_FE = 5;  % kernel extent in frequency-encoded direction
  Nsrc_PE = 5;  % kernel extent in phase-encoded direction (a.k.a., Nblocks)

  % an even number of kernels might make sense here. if there were two,
  % there could be one "looking forward" and one "looking back".

  Nkernel = 4;  % number of shifted kernels: 1 <= Nkernel < Nblocks

  
  %==--------------------------------------------------------------------==%

  Ncol = size(k_reduFOV, 1);
  Ncha = size(k_reduFOV, 3);

  k_fullFOV = zeros(Ncol, Nlin, Ncha);

  % collect lines starting at first fullFOV k-space line
  %  (1-based)
  index_datlines = [0:R:(Nlin-1)] + 1;
  token_datlines = zeros(1,Nlin);
  token_datlines(index_datlines) = 2;

  index_full2dat = zeros(1,Nlin);
  index_full2dat(index_datlines) = 1:length(index_datlines);


  k_fullFOV(:,index_datlines,:) = k_reduFOV;

  kernels = [ones(Nsrc_PE,1)  *     (0:Nsrc_PE-1)] ...
            - [(0:Nsrc_PE-1)' * ones(1,Nsrc_PE)];

  % a lazy way to re-order the kernels inside-out:
  [sorted, col_permute] = sort(abs( ([-(Nsrc_PE-1)/2]:[+(Nsrc_PE-1)/2]) - 0.25 ) );
  kernels = kernels(:,col_permute);
  
  width_FE = (Nsrc_FE-1)/2;
  width_PE = (Nsrc_PE-1)/2;

  kernel_ind_FE = (-width_FE : +width_FE);

  for ind = 1:Nkernel,

    if ( Nkernel == 1 ),
      kernel_ind_PE = (-width_PE : +width_PE) * R;
    else,
      kernel_ind_PE = kernels(:,ind).' * R;
    end;
    
    for ii = 1:Ncol,
      for jj = 1:(Nlin/R),

        Xpos = ii;

        Ypos = index_datlines(jj);

        Xsrc = mod([ (Xpos-1) + kernel_ind_FE ], Ncol) + 1;
        Ysrc = mod([ (Ypos-1) + kernel_ind_PE ], Nlin) + 1;

        Xtrg = ii;
        Ytrg = (Ypos+1):(Ypos+R-1);

        s_src = mrir_array_GRAPPA_1d__extract_entries(k_reduFOV, Xsrc, index_full2dat(Ysrc), 1:Ncha);
	s_trg = G{ind} * s_src;
	
	k_fullFOV = mrir_array_GRAPPA_1d__insert_entries(k_fullFOV, s_trg, Xtrg, Ytrg, 1:Ncha, Nkernel);

      end;
    end;
    
  end;

  
  %==--------------------------------------------------------------------==%

  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total read time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));

  
  if ( nargout > 0 ),
    varargout{1} = k_fullFOV;
  end;
  

  return;



%**************************************************************************%
function [m, X, Y, Z] = mrir_array_GRAPPA_1d__extract_entries(M, x, y, z)


  [X, Y, Z] = meshgrid(x,y,z);
  ind = sub2ind(size(M), X, Y, Z);
  m = reshape(M(ind), [], 1);


  return;



%**************************************************************************%
function M = mrir_array_GRAPPA_1d__insert_entries(M, v, x, y, c, N)


  [X, Y, C] = ndgrid(x,y,c);
  ind = sub2ind(size(M), X, Y, C);
  M(ind) = M(ind(:)) + (v/N);


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_1d_recon.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
