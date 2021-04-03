



      e = eig(errcovmtx_inv);
      emax = max(e);
      emin = min(e);
      
      f = 0.00005;

      lambda = (emax-emin)*f + emin;
      
      
      % kellman2001SENSE
      % sodickson2000tailored
      
%      errcovmtx = inv(errcovmtx_inv + 1e-8*eye(size(errcovmtx_inv)));
      errcovmtx = inv( errcovmtx_inv + lambda*eye(size(errcovmtx_inv)) );

      A = diag(1./diag(errcovmtx * E' * covmtxinv * E));

      U = A * errcovmtx * E' * covmtxinv;
      rho = U * E;
      
      nu = rho * errcovmtx * A;
      
      
      