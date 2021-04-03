function [M_phase_lumped, M_phase_coax] = mrir_array_modemix_bin_phz(M, phase_trace, phase_circuit)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/06
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
    
  % subtract board traces from matrix phases
  M_phase_shifter = rad2deg(angle(  M ./ exp(i*repmat(deg2rad(phase_board_trace), [1, 16]))  ));

  phase_circuit_mat = repmat(reshape(phase_circuit, 1, 1, []), [Nchan, Nmode, 1]);


  % first calculate discrepancy between the phase shifters we will add to the
  % board and the phases implementable by the lumped version of the phase
  % shifter
  phase_diff = phase_circuit_mat - repmat(modemixmatrix_phase_shifter, [1 1 length(phase_circuit)]);


  phase_diff(find(phase_diff < 0)) = Inf;

  % find the lumped phase shifter circuit that is closest (and above!) the desired phase shift
  [buffer, circuit_index] = min(phase_diff, [], 3);

  % set 233 degrees to the equivalent -127 degrees
  phase_circuit_wrap = [phase_circuit(1) - 360, phase_circuit(2:end)];

  % extract the lumped phase shifter circuit for each channel in the matrix
  M_phase_lumped = phase_circuit(circuit_index);


  % the phase shift implemented by coax is the remainder of the required phase shifter
  M_phase_coax = M_phase_lumped - M_phase_shifter;


  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
