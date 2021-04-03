function A = calmat(C, k, r, t)
  % k -> Spatial kernel dimension. [7 7 1]
  % r-> crop [360 360 1]
  % t -> Time kernel dimension.  [1]
  % C -> kspace 

  [sx, sy, sz, nc, ec, mm] = size(C);

  assert(numel(k) == numel(r));
  if (numel(k) == 1)
    r = [min(r,   sx), min(r,   sy), min(r,   sz)];
    k = [min(k, r(1)), min(k, r(2)), min(k, r(3))];
  end

  assert(mod(r(1), k(1)) == 0);
  assert(mod(r(2), k(2)) == 0);
  assert(mod(r(3), k(3)) == 0);

  A = single(zeros(max(1, ec - t + 1) * prod(r(r > 1) - k(r > 1) + 1), mm * t * nc * prod(k(k > 1))));

  ctr = 0;
  for ee = 1:max(min(ec - t + 1, ec), 1)
    for xx = 1:max(min(r(1) - k(1) + 1, r(1)), 1)
      for yy = 1:max(min(r(2) - k(2) + 1, r(2)), 1)
        for zz = 1:max(min(r(3) - k(3) + 1, r(3)), 1)

          erange = blockrange(ee, ec, t);
          xrange = blockrange(xx, r(1), k(1));
          yrange = blockrange(yy, r(2), k(2));
          zrange = blockrange(zz, r(3), k(3));

          block = C(xrange, yrange, zrange, :, erange, :);

          ctr = ctr + 1;

          A(ctr, :) = block(:);

        end
      end
    end
  end

end

function range = blockrange(startnum, bound, blocksize)
  range = (startnum - 1) + [1:blocksize];
  range = range(range <= bound);
end
