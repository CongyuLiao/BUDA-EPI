function [M_bin, tally] = mrir_array_modemix_bin_phz_inductors(M, phase_shift_deg)

  

  phase_shift = deg2rad(phase_shift_deg);
  M_shift = exp(-i * repmat(phase_shift, [1, size(M,2)]));
  
  M_phz = angle(M .* M_shift);
  
%  edges_deg = [-180:45:-45, -5,+5, +45:45:+180];
  edges_deg = [-180:45:+180];
  edges = deg2rad(edges_deg);

  [count, bin] = histc(M_phz, edges);

  phz_levels = edges + pi/4;
  phz_bin = phz_levels(bin);

  M_bin = abs(M) .* exp(i*phz_bin);

  
  
  for ind = 1:length(phz_levels)-1,
    
    tally(ind) = sum( phz_bin(:) == phz_levels(ind) );
    
  end;