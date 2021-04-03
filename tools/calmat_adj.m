function C = calmat_adj(A, k, r, t, calibsize)
% calibsize  = size of C

  sx = calibsize(1);
  sy = calibsize(2);
  sz = calibsize(3);
  nc = calibsize(4);
  ec = calibsize(5);
  mm = calibsize(6);

  C = single(zeros(sx, sy, sz, nc, ec, mm));
  D = single(zeros(sx, sy, sz, nc, ec, mm));
  O = single(ones(numel(blockrange(1, sx, k)), ...
                  numel(blockrange(1, sy, k)), ...
                  numel(blockrange(1, sz, k)), ...
                  nc,                          ...
                  numel(blockrange(1, ec, t)), ...
                  mm                         ));

  ctr = 0;
  for ee = 1:min(ec - t + 1, ec)
    for xx = 1:min(r - k + 1, sx)
      for yy = 1:min(r - k + 1, sy)
        for zz = 1:min(r - k + 1, sz)

          ctr = ctr + 1;

          erange = blockrange(ee, ec, t);
          xrange = blockrange(xx, sx, k);
          yrange = blockrange(yy, sy, k);
          zrange = blockrange(zz, sz, k);

          block = reshape(A(ctr, :), numel(xrange), numel(yrange), numel(zrange), nc, numel(erange), mm);

          C(xrange, yrange, zrange, :, erange, :) = C(xrange, yrange, zrange, :, erange, :) + block;
          D(xrange, yrange, zrange, :, erange, :) = D(xrange, yrange, zrange, :, erange, :) + O;

        end
      end
    end
  end

  C = C./D;

end

function range = blockrange(startnum, bound, blocksize)
  range = (startnum - 1) + [1:blocksize];
  range = range(range <= bound);
end
