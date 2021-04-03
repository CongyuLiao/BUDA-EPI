function varargout = mrir_array_GRAPPA_1d(k_reduFOV, k_ACS, R, varargin)
%MRIR_ARRAY_GRAPPA_1D
%
% G = mrir_array_GRAPPA_1d(k_reduFOV, k_ACS, R)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/10
% $Id: mrir_array_GRAPPA_1d_kernel.m,v 1.3 2008/04/01 05:53:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG; if ( isempty(DEBUG) ), DEBUG = 0; end;


  %==--------------------------------------------------------------------==%

  t0 = clock;

  
  if ( nargin >= 4 ),
    index_datlines = varargin{4-3};
  else,
    index_datlines = [];
  end;

  if ( nargin >= 5 ),
    index_acslines = varargin{5-3};
  else,
    index_acslines = [];
  end;


  %%% kernel parameters:

  if ( nargin >= 6 ),
    Nsrc_FE  = varargin{6-3};
  else,
    Nsrc_FE = 5;    % kernel extent in frequency-encoded direction
  end;

  if ( nargin >= 7 ),
    Nsrc_PE  = varargin{7-3};
  else,
    Nsrc_PE = 5;    % kernel extent in phase-encoded direction (a.k.a., Nblocks)
  end;


  % an even number of kernels might make sense here. if there were two,
  % there could be one "looking forward" and one "looking back".

  if ( nargin >= 8 ),
    Nkernel  = varargin{8-3};
  else,
    Nkernel  = 4;   % number of shifted kernels: 1 <= Nkernel < Nblocks
  end;


  %==--------------------------------------------------------------------==%

  % decisions:

  % since ACS lines usually overlap with collected lines, should the
  % kernel be fit to only ACS lines or should the data lines be
  % substituted? according to breuer2004autocalibrated, the ACS lines do
  % not need to have the same contrast, so probably weights should be
  % computed only from ACS lines.

  % but if kernel is large, will extend way beyond ACS lines. conclusion:
  % either use only ACS lines or only data lines for sources.


  Ncol = size(k_reduFOV, 1);
  Nlin = size(k_reduFOV, 2) * R;
  Ncha = size(k_reduFOV, 3);

  % number of collected ACS lines
  Nacs = size(k_ACS, 2);

  % number of fits of the GRAPPA kernel obtained from this data
  Nfit = Ncol * Nacs / R;


  %==--------------------------------------------------------------------==%

  % if we are not only using ACS lines to compute the weights, then we
  % define reference lines as those that are not collected along with the
  % normal data, i.e., acslines but not datlines.

  if ( isempty(index_datlines) ),
    % collect lines starting at first fullFOV k-space line
    %  (1-based)
    index_datlines = [0:R:(Nlin-1)] + 1;
  end;
  token_datlines = zeros(1,Nlin);
  token_datlines(index_datlines) = 2;

  index_full2dat = zeros(1,Nlin);
  index_full2dat(index_datlines) = 1:length(index_datlines);


  if ( isempty(index_acslines) ),
    %  (1-based)
    index_acslines = [ (Nlin/2 - Nacs/2) : (Nlin/2 + Nacs/2)-1 ] + 1;
  end;
  token_acslines = zeros(1,Nlin);
  token_acslines(index_acslines) = 1;

  index_full2acs = zeros(1,Nlin);
  index_full2acs(index_acslines) = 1:length(index_acslines);

  token_xxxlines = token_datlines .* token_acslines;
  token_reflines = token_acslines .* (~token_datlines);
  %  (1-based)
  index_reflines = find(token_reflines);

  % number of true reference lines
  Nref = size(index_reflines, 2);


  %==--------------------------------------------------------------------==%

  kernels = [ones(Nsrc_PE,1)  *     (0:Nsrc_PE-1)] ...
            - [(0:Nsrc_PE-1)' * ones(1,Nsrc_PE)];

  % a lazy way to re-order the kernels inside-out:
  [sorted, col_permute] = sort(abs( ([-(Nsrc_PE-1)/2]:[+(Nsrc_PE-1)/2]) - 0.25 ) );
  kernels = kernels(:,col_permute);

  width_FE = (Nsrc_FE-1)/2;
  width_PE = (Nsrc_PE-1)/2;

  ind_FE = mod([ -width_FE : (Ncol-1) + width_FE ], Ncol) + 1;
  kernel_ind_FE = (-width_FE : +width_FE);
    
  for ind = 1:Nkernel,
    
    if ( Nkernel == 1 ),
      kernel_ind_PE = (-width_PE : +width_PE) * R;
    else,
      kernel_ind_PE = kernels(:,ind).' * R;
    end;

    % preallocate for time efficiency
    S_src = zeros(Ncha * Nsrc_FE * Nsrc_PE, Nfit);
    S_trg = zeros(Ncha * (R-1),             Nfit);

    for ii = 1:Ncol,
      for jj = 1:(Nacs/R),

        Xpos = ii;

        % find data line immediately below (for default kernel) current
        % ACS line
        Ypos = index_datlines(max(find((index_reflines(jj*(R-1)) - index_datlines) > 0)));

        Xsrc = mod([ (Xpos-1) + kernel_ind_FE ], Ncol) + 1;
        Ysrc = mod([ (Ypos-1) + kernel_ind_PE ], Nlin) + 1;

        Xtrg = ii;
        Ytrg = (Ypos+1):(Ypos+R-1);

        [s_src, X1,Y1] = mrir_array_GRAPPA__extract_entries(k_reduFOV, Xsrc, index_full2dat(Ysrc), 1:Ncha);
        [s_trg, X2,Y2] = mrir_array_GRAPPA__extract_entries(k_ACS,     Xtrg, index_full2acs(Ytrg), 1:Ncha);

        if ( DEBUG ),
          mrir_array_GRAPPA_1d__plot_kernel(X1, Y1, X2, Y2, Ncol, Nlin, ...
                                            Nsrc_FE, Nsrc_PE, ii, jj);
        end;

        % S_src: [Ncha * Nsrc_FE * Nsrc_PE] x Nfit
        S_src(:, sub2ind([Nacs/R, Ncol], jj, ii)) = s_src;
	
        % S_trg: [Ncha * (R-1)] x Nfit
        S_trg(:, sub2ind([Nacs/R, Ncol], jj, ii)) = s_trg;

      end;
    end;

    % G:  [Ncha * (R-1)] x [Ncha * Nsrc_FE * Nsrc_PE]
    G{ind} = S_trg * pinv(S_src);

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
    varargout{1} = G;
  end;


  return;



%**************************************************************************%
function [m, X, Y, C] = mrir_array_GRAPPA__extract_entries(M, x, y, c)


  [X, Y, C] = ndgrid(x,y,c);
  ind = sub2ind(size(M), X, Y, C);
  m = reshape(M(ind), [], 1);


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_1d_kernel.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
