function [kernel_conv, kernel_corr] = mrir_array_GRAPPA_conv_kernel(G)
%MRIR_ARRAY_GRAPPA_CONV_KERNEL
%
% kernel = mrir_array_GRAPPA_conv_kernel(G)
%
%
%   k:  [Nsrcx, Nsrcy, Nsrcz, Nsrcc, Ntrgc]
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/05
% $Id: mrir_array_GRAPPA_conv_kernel.m,v 1.4 2008/12/14 04:30:07 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  kernels1 = mrir_array_GRAPPA__pick_kernels(G.Nsrcy);
  kernels2 = mrir_array_GRAPPA__pick_kernels(G.Nsrcz);

  kernel_y = kernels1(:,1).' * G.R1;
  kernel_z = kernels2(:,1).' * G.R2;


  % G:  [NCha * (R1-1) * (R2-1)] x [NCha * Nsrcx * Nsrcy * Nsrcz]



  Nx = G.Nsrcx;
  Ny = max([ 1,   max(kernel_y) - min(kernel_y) + (G.R1-1)   ]);
  Nz = max([ 1,   max(kernel_z) - min(kernel_z) + (G.R2-1)   ]);

  NCha = size(G.kernel{1}, 2) / G.Nsrcx / G.Nsrcy / G.Nsrcz;
  
  R1 = G.R1;
  R2 = G.R2;

  Ry = max([1, (R1-1)]);
  Rz = max([1, (R2-1)]);
  
  
  kernel_corr = zeros(Nx, Ny, Nz, NCha, NCha);

  
  % assume Nx is always odd-valued
  Xcenter = ceil(Nx/2);

  Ypos = 1 - min(kernel_y);
  Ycenter = Ypos + (R1-1);

  Zpos = 1 - min(kernel_z);
  Zcenter = Zpos + (R2-1);

  for cha = 1:NCha,
    kernel_corr(Xcenter, Ycenter, Zcenter, cha, cha) = 1;
  end;
  
  ind = 0;

  for ind_z = 1:Rz,
    for ind_y = 1:Ry,

      ind = ind + 1;
      
      % extract all channels for this target point, w: NCha x [NCha * Nsrcx * Nsrcy * Nsrcz]
      w = G.kernel{1}(ind:(Ry*Rz):end,:);
      
      W = reshape( w.', G.Nsrcx, G.Nsrcy, G.Nsrcz, NCha, NCha);

      %%%% test:
      %% W1 = reshape(w(1,:), G.Nsrcx, G.Nsrcy, G.Nsrcz, NCha);
      %% W(:,:,:,:,1) - W1

      
      
      kernel_corr(:, kernel_y + Ypos + Ry - ind_y, kernel_z + Zpos + Rz - ind_z, :, :) = W;

    end;
  end;

  % definition of convolution includes reversing the correlation kernel along each convolution direction
  kernel_conv = flipdim(flipdim(flipdim(kernel_corr, mrir_DIM_COL), mrir_DIM_LIN), mrir_DIM_PAR);
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_GRAPPA_conv_kernel.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
