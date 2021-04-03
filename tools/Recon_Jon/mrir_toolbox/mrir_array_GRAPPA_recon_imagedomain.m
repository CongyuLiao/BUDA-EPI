function [i_conv, k_conv, K] = mrir_array_GRAPPA_recon_imagedomain(k_redu, K, evp, varargin)
%MRIR_ARRAY_GRAPPA_RECON_IMAGEDOMAIN
%
% [i_conv, k_conv, K] = mrir_array_GRAPPA_recon_imagedomain(raw, G, evp)
%
% NOTE: raw represents full k-space data (containing zeroed-out lines)
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/15
% $Id: mrir_array_GRAPPA_recon_imagedomain.m,v 1.1 2008/12/16 00:24:53 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  t0 = clock;


  %==--------------------------------------------------------------------==%

  NCol = mrir_ice_dimensions(k_redu, 'col');
  NLin = mrir_ice_dimensions(k_redu, 'lin');
  NCha = mrir_ice_dimensions(k_redu, 'cha');
  NPar = mrir_ice_dimensions(k_redu, 'par');

  flReadoutOSFactor = 2;

  
  % its AMAZING that this works...
  k_datlines1 = evp.NFirstLin : K.R1 : evp.NLinMeas;
  full2dat1(k_datlines1) = 1:length(k_datlines1);
  if ( length(k_datlines1) ~= evp.RawLin ), warning('i am confused! (1)'); end;

  k_datlines2 = evp.NFirstPar : K.R2 : evp.NParMeas;
  full2dat2(k_datlines2) = 1:length(k_datlines2);
  if ( length(k_datlines2) ~= evp.RawPar ), warning('i am confused! (2)'); end;

  fprintf(1, '<i> [%s]: allocating...', mfilename);
  k_fill = complex(zeros( NCol, evp.NLinMeas, NCha, 1, 1, 1, 1, 1, evp.NParMeas, class(k_redu) ));
    
  fprintf(1, 'and...');
  k_fill(:, k_datlines1, :, 1,1,1,1,1, k_datlines2) = k_redu;
  fprintf(1, 'done!\n');

    
  %i_redu = mrir_conventional(k_fill, 0);
% Kawin
  i_redu = mrir_iDFT(mrir_iDFT(k_fill,1),2);

  %==--------------------------------------------------------------------==%

  t1 = clock;

  if ( isstruct(K) ),
    G = K;
    k = mrir_array_GRAPPA_conv_kernel(G);
    K = mrir_iDFT(mrir_iDFT(mrir_iDFT(k, mrir_DIM_COL, NCol), mrir_DIM_LIN, evp.NLinMeas), mrir_DIM_PAR, evp.NParMeas);
  end;


  for trg = 1:NCha,
    i_conv(:,:,trg) = sum( [i_redu .* (squeeze(K(:,:,1,:,trg)))], mrir_DIM_CHA);
  end

  t2 = clock;
  
  
  k_conv = mrir_fDFT_freqencode(mrir_fDFT_phasencode(i_conv));
  
  
  runtime_seconds = etime(t2,t1);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total recon time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  %==--------------------------------------------------------------------==%

  tN = clock;
  runtime_seconds = etime(tN,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total execution time = %s', TIME);
%  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_recon_imagedomain.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
