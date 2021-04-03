function kernels = mrir_array_GRAPPA__pick_kernels(Nsrc)


  kernels = [ones(Nsrc,1)  *     (0:Nsrc-1)] ...
            - [(0:Nsrc-1)' * ones(1,Nsrc)];

  % a lazy way to re-order the kernels inside-out:
  [sorted, col_permute] = sort(abs( ([-(Nsrc-1)/2]:[+(Nsrc-1)/2]) - 0.25 ) );
  kernels = sortrows(kernels(:,col_permute));


  return;