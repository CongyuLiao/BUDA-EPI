function fig = mrir_array_GRAPPA_plottest(NFirstLin, NLinMeas, NImageLins, ...
					  NAFLin)


  
  fig = figure('name', mfilename); axis; hold on;
  
  
  datlines1 = NFirstLin : NAFLin : NLinMeas;

  acclines1 = setdiff(1:NLinMeas, datlines1);
  
  x = [1,256];
  
  
  for ind = 1:length(datlines1),   
    plot(x, repmat(datlines1(ind), 2), 'k-');
  end;
  
  for ind = 1:length(acclines1),   
    plot(x, repmat(acclines1(ind), 2), 'k--');
  end;
  
  
  axis([-10 266 -5 133]);
  axis equal;
  set(gca, 'XTick', [0,64:64:256])
  set(gca, 'YTick', [0,32:32:128])
