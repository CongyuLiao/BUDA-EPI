function mrir_array_combine_noise_correction_plots(SNRmax, Ncha)



  snr = linspace(0, SNRmax, 101);

  for ind = 1:length(Ncha),

    [Am_per_Sm{ind}, At_per_St{ind}, Sm_per_St{ind}, Am_per_St{ind}] = ...
        mrir_array_combine_noise_correction_lookup(snr, Ncha(ind));

    Rn{ind} = At_per_St{ind} ./ Am_per_Sm{ind};
    Dn{ind} = Am_per_Sm{ind}  - At_per_St{ind};

  end;

  save(mfilename)
  
  for ind = length(Ncha): -1 : 1,
  
    if ( isempty(Am_per_Sm{ind}) ),
      Am_per_Sm = {Am_per_Sm{1:ind-1}};
      At_per_St = {At_per_St{1:ind-1}};
      Sm_per_St = {Sm_per_St{1:ind-1}};
      Am_per_St = {Am_per_St{1:ind-1}};
      Ncha = Ncha(1:ind-1);
    end;
    
  end;
    

  cmap = jet(length(Ncha));
    
  
  %==-----------------------------------------------------------------------==%

  fig_a = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(At_per_St{ind}, Am_per_St{ind}, 'color', cmap(ind, :));
  end;
  plot([0, SNRmax], [0, SNRmax], 'k--');
  axis equal;
  set(gca, 'XLim', [0, SNRmax], 'YLim', [0, SNRmax]);
  xlabel('A_n / \sigma');
  ylabel('M_n / \sigma');
  title('constantinides, fig 2(a)');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'SouthEast');
  
  fig_b = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(At_per_St{ind}, Sm_per_St{ind}, 'color', cmap(ind, :));
  end;
  plot([0, 25], [1, 1], 'k--');
  set(gca, 'YLim', [0.50, 1.10]);
  xlabel('A_n / \sigma');
  ylabel('\sigma_M_n / \sigma');
  title('constantinides, fig 2(a)');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'SouthEast');

  fig_c = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(Am_per_St{ind}, Am_per_St{ind} - At_per_St{ind}, 'color', cmap(ind, :));
  end;
  plot([0, 25], [0, 25], 'k--');
  axis equal;
  xlabel('"measured" SNR (M_n / \sigma)');
  ylabel('subtractive correction factor ((M_n - A_n) / \sigma)');
  title('constantinides, fig 3');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'NorthEast');


  %==-----------------------------------------------------------------------==%

  fig_d = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(At_per_St{ind}, Am_per_Sm{ind}, 'color', cmap(ind, :));
  end;
  plot([0, SNRmax], [0, SNRmax], 'k--');
  axis equal;
  set(gca, 'XLim', [0, SNRmax], 'YLim', [0, SNRmax]);
  xlabel('"true" SNR');
  ylabel('"measured" SNR');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'SouthEast');
  box on;
  title('bias in SNR measured from magnitude images');


  fig_e = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(Am_per_Sm{ind}, Rn{ind}, 'color', cmap(ind, :));
  end;
  plot([0, 1.1*SNRmax], [1, 1], 'k--');
  set(gca, 'XLim', [0, SNRmax], 'YLim', [0, 1.1]);
  xlabel('"measured" SNR');
  ylabel('divisive correction factor (true / meas)');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'SouthEast');
  box on;
  title('SNR bias correction from non-central \chi distribution');



  fig_f = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph(ind) = plot(Am_per_Sm{ind}, Dn{ind}, 'color', cmap(ind, :));
  end;
  plot([0, SNRmax], [0, SNRmax], 'k--');
  axis equal;
  set(gca, 'XLim', [0, SNRmax], 'YLim', [0, SNRmax]);
  xlabel('"measured" SNR');
  ylabel('subtractive correction factor (meas - true)');
  legend(fliplr(ph), flipud(cellstr(num2str(Ncha(:)))), 'Location', 'NorthWest');
  box on;
  title('SNR bias correction from non-central \chi distribution');

  
  
  fig_g = figure('name', mfilename); axis; hold on;
  for ind = 1:length(Ncha),
    ph1 = plot(Am_per_Sm{ind}, Dn{ind}, 'k');
    ph2 = plot(Am_per_St{ind}, Am_per_St{ind} - At_per_St{ind}, 'r');
  end;
  plot([0, SNRmax], [0, SNRmax], 'k--');
  axis equal;
  set(gca, 'XLim', [0, SNRmax], 'YLim', [0, SNRmax]);
  xlabel('"measured" SNR');
  ylabel('subtractive correction factor');
  legend([ph1, ph2], 'original', 'corrected', 'Location', 'NorthWest');
  box on;
  title('SNR correction differs from "mean correction in \sigma units"');

  