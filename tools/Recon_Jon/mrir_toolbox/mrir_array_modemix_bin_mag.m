function [A_attenuation_log_dB, M_attenuation_log] = mrir_array_modemix_bin_mag(M, FLAG__Normalize_Individual_Modes)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/oct/25
% $Id: mrir_array_modemix_bin_mag.m,v 1.1 2008/11/06 05:45:06 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  step_attenuator_chip = 15;

  
  if ( step_attenuator_chip == 15 ),

    % DAT-15R5-PN: 0.5 to 15.5 dB of *power* attenuation

    attenuation_steps = sort([0, -0.5, -1, -2, -4, -8]);
    attenuation_levels_dB = sum(attenuation_steps):0.5:attenuation_steps(end);

  end;

  if ( step_attenuator_chip == 31 ),

    % DAT-31-PN: 1 to 31 dB of *power* attenuation

    attenuation_31dB_steps = sort([-1, -2, -4, -8, -16]);
    attenuation_levels_dB = sum(attenuation_31dB_steps):1.0:attenuation_31dB_steps(end);

  end;

  attenuation_levels = 10.^(attenuation_levels_dB/20);

    
  A = abs(M);

  if ( FLAG__Normalize_Individual_Modes ),

    amin = min(A, [], 1);
    Amin = repmat(amin, size(A, 1), 1);

    amax = max(A, [], 1);
    Amax = repmat(amax, size(A, 1), 1);
    
    A = ( A - Amin ) ./ ( Amax - Amin );

    attenuation_levels_linear = 10.^( attenuation_levels_dB/20 );
    
    A = [A * ( attenuation_levels_linear(end) - attenuation_levels_linear(1) )] + attenuation_levels_linear(1);
    
    
  else,

    Amin = min(A(:));
    Amax = max(A(:));
    
    % is it OK to subtract this constant from the matrix? (is constant vector in nullspace?) ((sleepy))
    A = ( A - Amin ) ./ ( Amax - Amin );
  
  end;
  

  A_dB = 20*log10(A);


  %---------------------------------------------------------------------------%

  edges_log = [-Inf attenuation_levels_dB + [diff(attenuation_levels_dB)/2, +Inf]];
  [count_log, bin_log] = histc(A_dB, edges_log);

  A_attenuation_log_dB = attenuation_levels_dB(bin_log);

  M_attenuation_log = ( 10.^(A_attenuation_log_dB/20) ) .* exp(i*angle(M));



  %---------------------------------------------------------------------------%

  edges_lin = [-Inf attenuation_levels + [diff(attenuation_levels)/2, +Inf]];
  [count_lin, bin_lin] = histc(A, edges_lin);

  A_attenuation_lin = attenuation_levels(bin_lin);

  M_attenuation_lin = A_attenuation_lin .* exp(i*angle(M));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_modemix_bin_mag.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
