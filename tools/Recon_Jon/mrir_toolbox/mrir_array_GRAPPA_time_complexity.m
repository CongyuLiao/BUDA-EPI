function [performance_elapsed, performance_loadavg] = mrir_array_GRAPPA_time_complexity(noisecov, NCol, NLin)
 
  NCha = unique(size(noisecov));
  
  k_rand = mrir_noise_synthesize([NCol, NLin, NCha], noisecov);
  
  
  NRefLin = 32;
  
  Nsrcx = 3;
  Nsrcy = 4;
  
  performance_elapsed = [];
  performance_loadavg = mrir_sysutil__loadavg;
  
  for R = 2:6,
   
  
  evp.NAFLin = R;
  evp.NFirstLin = 1;
  evp.NLinMeas = NLin;
  
  evp.RawLin = ceil(NLin / R);

  evp.NFirstRefLin = (NLin/2) - (NRefLin/2);
  evp.NRefLin = NRefLin;

  evp
  
  
  [k_redu, k_ACS] = mrir_array_accelerated_synthesize(k_rand, evp, 1);

    
    
    for cha = 4:96,

      t1 = clock;
      
      G = mrir_array_GRAPPA_2d_kernel(k_ACS(:,:,1:cha), R, 1, Nsrcx, Nsrcy, 1, 1, 1);

      t2 = clock;

      
      dat_R4_full = mrir_array_GRAPPA_2d_recon(k_redu(:,:,1:cha), G, ...
                                               R, ...
                                               1, ...
                                               NLin, ...
                                               1, ...
                                               1, ...
                                               1);

      t3 = clock;
      
      img_R4_full = mrir_conventional(dat_R4_full);
      rss_R4_full = mrir_array_combine_rss(img_R4_full);

      t4 = clock;
      
      performance_elapsed = cat(1, performance_elapsed, [etime(t2, t1), etime(t3, t2), etime(t4, t3), etime(t4, t1)]);
      performance_loadavg = cat(1, performance_loadavg, mrir_sysutil__loadavg);
      
    end;
    
  
    
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
  
  
  
%    i_rand = mrir_conventional(k_rand, struct('flReadoutOSFactor', 1, 'lPhaseEncodingLines', NLin));

  